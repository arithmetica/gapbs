
/* Taken from Galois 
  Nov 28 2018
 */

#include <cassert>

#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

#include "file_graph.h"

FileGraph::FileGraph()
  : masterMapping(0), masterLength(0), masterFD(0),
    outIdx(0), outs(0), edgeData(0),
    numEdges(0), numNodes(0)
{
}

FileGraph::~FileGraph() {
  if (masterMapping)
    munmap(masterMapping, masterLength);
  if (masterFD)
    close(masterFD);
}

void FileGraph::parse(void* m) {
  //parse file
  uint64_t* fptr = (uint64_t*)m;
  uint64_t version = convert_le64(*fptr++);
  if (version != 1)
    std::cout << "unknown file version: " << version << "\n";
  sizeofEdge = convert_le64(*fptr++);
  numNodes = convert_le64(*fptr++);
  numEdges = convert_le64(*fptr++);
  outIdx = fptr;
  fptr += numNodes;
  uint32_t* fptr32 = (uint32_t*)fptr;
  outs = fptr32; 
  fptr32 += numEdges;
  if (numEdges % 2)
    fptr32 += 1;
  edgeData = (char*)fptr32;
}


void FileGraph::structureFromMem(void* mem, size_t len, bool clone) {
  masterLength = len;

  if (clone) {
    int _MAP_BASE = MAP_ANONYMOUS | MAP_PRIVATE;
#ifdef MAP_POPULATE
    _MAP_BASE |= MAP_POPULATE;
#endif
    
    void* m = mmap(0, masterLength, PROT_READ | PROT_WRITE, _MAP_BASE, -1, 0);
    if (m == MAP_FAILED) {
      std::cout << "failed copying graph" << "\n";
    }
    memcpy(m, mem, len);
    parse(m);
    masterMapping = m;
  } else {
    parse(mem);
    masterMapping = mem;
  }
}

void* FileGraph::structureFromGraph(FileGraph& g, size_t sizeof_edge_data) {
  // Allocate
  size_t common = g.masterLength - (g.sizeofEdge * g.numEdges);
  size_t len = common + (sizeof_edge_data * g.numEdges);
  int _MAP_BASE = MAP_ANONYMOUS | MAP_PRIVATE;
#ifdef MAP_POPULATE
  _MAP_BASE |= MAP_POPULATE;
#endif
  void* m = mmap(0, len, PROT_READ | PROT_WRITE, _MAP_BASE, -1, 0);
  if (m == MAP_FAILED) {
    std::cout << "failed copying graph" << "\n";
  }
  memcpy(m, g.masterMapping, common);
  uint64_t* fptr = (uint64_t*)m;
  fptr[1] = convert_le64(sizeof_edge_data);
  structureFromMem(m, len, false);

  return edgeData;
}

void* FileGraph::structureFromArrays(uint64_t* out_idx, uint64_t num_nodes,
      uint32_t* outs, uint64_t num_edges, size_t sizeof_edge_data) {
  uint64_t nBytes = sizeof(uint64_t) * 4; // version, sizeof_edge_data, numNodes, numEdges

  nBytes += sizeof(uint64_t) * num_nodes;
  nBytes += sizeof(uint32_t) * num_edges;
  if (num_edges % 2)
    nBytes += sizeof(uint32_t); // padding
  nBytes += sizeof_edge_data * num_edges;
 
  int _MAP_BASE = MAP_ANONYMOUS | MAP_PRIVATE;
#ifdef MAP_POPULATE
  _MAP_BASE |= MAP_POPULATE;
#endif
  
  char* base = (char*) mmap(0, nBytes, PROT_READ | PROT_WRITE, _MAP_BASE, -1, 0);
  if (base == MAP_FAILED) {
    base = 0;
    std::cout << "failed allocating graph" << "\n";
  }
  
  uint64_t* fptr = (uint64_t*) base;
  *fptr++ = convert_le64(1);
  *fptr++ = convert_le64(sizeof_edge_data);
  *fptr++ = convert_le64(num_nodes);
  *fptr++ = convert_le64(num_edges);

  for (size_t i = 0; i < num_nodes; ++i)
    *fptr++ = convert_le64(out_idx[i]);
  uint32_t* fptr32 = (uint32_t*) fptr;
  for (size_t i = 0; i < num_edges; ++i)
    *fptr32++ = convert_le32(outs[i]);

  structureFromMem(base, nBytes, false);
  return edgeData;
}

void FileGraph::structureFromFile(const std::string& filename, bool preFault) {
  masterFD = open(filename.c_str(), O_RDONLY);
  if (masterFD == -1) {
    std::cout << "failed opening " << filename << "\n";
  }

  struct stat buf;
  int f = fstat(masterFD, &buf);
  if (f == -1) {
    std::cout << "failed reading " << filename << "\n";
  }
  masterLength = buf.st_size;

  int _MAP_BASE = MAP_PRIVATE;
#ifdef MAP_POPULATE
  if (preFault)
    _MAP_BASE |= MAP_POPULATE;
#endif
  
  void* m = mmap(0, masterLength, PROT_READ, _MAP_BASE, masterFD, 0);
  if (m == MAP_FAILED) {
    m = 0;
    std::cout << "failed reading " << filename << "\n";
  }
  parse(m);
  masterMapping = m;

#ifndef MAP_POPULATE
  if (preFault) {
    Runtime::MM::pageIn(m, masterLength);
  }
#endif
}


void FileGraph::swap(FileGraph& other) {
  std::swap(masterMapping, other.masterMapping);
  std::swap(masterLength, other.masterLength);
  std::swap(sizeofEdge, other.sizeofEdge);
  std::swap(masterFD, other.masterFD);
  std::swap(outIdx, other.outIdx);
  std::swap(outs, other.outs);
  std::swap(edgeData, other.edgeData);
  std::swap(numEdges, other.numEdges);
  std::swap(numNodes, other.numNodes);
}

void FileGraph::cloneFrom(FileGraph& other) {
  structureFromMem(other.masterMapping, other.masterLength, true);
}

uint64_t FileGraph::getEdgeIdx(GraphNode src, GraphNode dst) const {
  for (uint32_t* ii = raw_neighbor_begin(src),
   *ee = raw_neighbor_end(src); ii != ee; ++ii)
    if (convert_le32(*ii) == dst)
      return std::distance(outs, ii);
  return ~static_cast<uint64_t>(0);
}

uint32_t* FileGraph::raw_neighbor_begin(GraphNode N) const {
  return (N == 0) ? &outs[0] : &outs[convert_le64(outIdx[N-1])];
}

uint32_t* FileGraph::raw_neighbor_end(GraphNode N) const {
  return &outs[convert_le64(outIdx[N])];
}

FileGraph::edge_iterator FileGraph::edge_begin(GraphNode N) const {
  return edge_iterator(N == 0 ? 0 : convert_le64(outIdx[N-1]));
}

FileGraph::edge_iterator FileGraph::edge_end(GraphNode N) const {
  return edge_iterator(convert_le64(outIdx[N]));
}

FileGraph::GraphNode FileGraph::getEdgeDst(edge_iterator it) const {
  return convert_le32(outs[*it]);
}

FileGraph::node_id_iterator FileGraph::node_id_begin() const {
  return boost::make_transform_iterator(&outs[0], Convert32());
}

FileGraph::node_id_iterator FileGraph::node_id_end() const {
  return boost::make_transform_iterator(&outs[numEdges], Convert32());
}

FileGraph::edge_id_iterator FileGraph::edge_id_begin() const {
  return boost::make_transform_iterator(&outIdx[0], Convert64());
}

FileGraph::edge_id_iterator FileGraph::edge_id_end() const {
  return boost::make_transform_iterator(&outIdx[numNodes], Convert64());
}

bool FileGraph::hasNeighbor(GraphNode N1, GraphNode N2) const {
  return getEdgeIdx(N1,N2) != ~static_cast<uint64_t>(0);
}

FileGraph::iterator FileGraph::begin() const {
  return iterator(0);
}

FileGraph::iterator FileGraph::end() const {
  return iterator(numNodes);
}

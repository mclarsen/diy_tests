#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include <diy/partners/swap.hpp>
#include <diy/reduce.hpp>
#include <diy/reduce-operations.hpp>

struct IntBlock
{
  int &m_int;
  IntBlock(int &an_int)
    : m_int(an_int)
  {
  }
};

struct AddIntBlock
{
  int &m_int;
  const diy::Master &m_master;

  AddIntBlock(diy::Master &master, int &an_int)
    : m_int(an_int),
      m_master(master)
  {
  }
  template<typename BoundsType, typename LinkType>                 
  void operator()(int gid,
                  const BoundsType &,  // local_bounds
                  const BoundsType &,  // local_with_ghost_bounds
                  const BoundsType &,  // domain_bounds
                  const LinkType &link) const
  {
    IntBlock *block = new IntBlock(m_int);
    LinkType *linked = new LinkType(link);
    diy::Master& master = const_cast<diy::Master&>(m_master);
    master.add(gid, block, linked);
  }
}; 

struct Redistribute
{
  typedef diy::RegularDecomposer<diy::DiscreteBounds> Decomposer;
  const diy::RegularDecomposer<diy::DiscreteBounds> &m_decomposer;
  Redistribute(const Decomposer &decomposer)
    : m_decomposer(decomposer)
  {}

  void operator()(void *v_block, const diy::ReduceProxy &proxy) const
  {
    IntBlock *block = static_cast<IntBlock*>(v_block);
    //
    // first round we have no incoming. Take the image we have,
    // chop it up into pieces, and send it to the domain resposible
    // for that portion 
    //
    const int world_size = m_decomposer.nblocks;

    if(proxy.in_link().size() == 0)
    {
      std::map<diy::BlockID, int> outgoing;

      for(int i = 0; i < world_size; ++i)
      {

        diy::BlockID dest = proxy.out_link().target(i); 

        outgoing[dest] = block->m_int; 
      } //for
      typename std::map<diy::BlockID,int>::iterator it;

      for(it = outgoing.begin(); it != outgoing.end(); ++it)
      {
        proxy.enqueue(it->first, it->second);
      }
    } // if
    else
    {
      // blend images according to vis order
      std::vector<int> ints;
      for(int i = 0; i < proxy.in_link().size(); ++i)
      {

        int incoming;
        int gid = proxy.in_link().target(i).gid;
        proxy.dequeue(gid, incoming); 
        ints.push_back(incoming);
      } // for

      for(int i = 0; i < (int)ints.size(); ++i)
      {
        std::cout<<ints[i] << " ";
      }
    } 

  } // operator
};

int main(void)
{
    MPI_Init(NULL, NULL);
    diy::mpi::communicator diy_comm;
    diy_comm = diy::mpi::communicator(MPI_COMM_WORLD);
    diy::DiscreteBounds global_bounds; 
 
    global_bounds.min[0] = 1;
    global_bounds.min[1] = 1;
    global_bounds.min[2] = 0;
    global_bounds.max[0] = 1024;
    global_bounds.max[1] = 1024;
    global_bounds.max[2] = 0;
    // tells diy to use all availible threads
    //const int num_threads = -1; 
    const int num_threads = 1; 
    const int num_blocks = diy_comm.size(); 
    const int magic_k = 8;

    diy::Master master(diy_comm, num_threads);

    diy::ContiguousAssigner assigner(num_blocks, num_blocks); 

    int data = diy_comm.rank();;
    AddIntBlock create(master, data);

    const int dims = 2;

    diy::RegularDecomposer<diy::DiscreteBounds> decomposer(dims, global_bounds, num_blocks);
    decomposer.decompose(diy_comm.rank(), assigner, create);

    diy::all_to_all(master, 
                    assigner, 
                    Redistribute(decomposer), 
                    magic_k);
}

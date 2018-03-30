#include <iostream>
#include <iterator>

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

void debug(void *b, 
           const diy::ReduceProxy &proxy,
           const diy::RegularSwapPartners &partners) 
{
  IntBlock *block = reinterpret_cast<IntBlock*>(b);
  unsigned int round = proxy.round();
  // count the number of incoming pixels
  if(proxy.in_link().size() > 0)
  {
      for(int i = 0; i < proxy.in_link().size(); ++i)
      {
        int gid = proxy.in_link().target(i).gid;
        if(gid == proxy.gid())
        {
          //skip revieving from self since we sent nothing
          continue;
        }
        int in;
        proxy.dequeue(gid, in);
      } // for in links
  } 

  if(proxy.out_link().size() == 0)
  {
    return;
  }
  // do compositing?? intermediate stage?
  const int group_size = proxy.out_link().size(); 
  const int current_dim = partners.dim(round);
  (void)current_dim; 
  std::vector<int> out(group_size);
  for(int i = 0; i < group_size; ++i)
  {
    out[i] = i;  
  } //for

  for(int i = 0; i < group_size; ++i)
  {
      if(proxy.out_link().target(i).gid == proxy.gid())
      {

      }
      else
      {
        proxy.enqueue(proxy.out_link().target(i), out[i]);
      }
  } //for 

} // reduce images

int main(void)
{
    MPI_Init(NULL, NULL);
    diy::mpi::communicator diy_comm;
    diy_comm = diy::mpi::communicator(MPI_COMM_WORLD);
    diy::DiscreteBounds global_bounds; 
 
    global_bounds.min[0] = 1;
    global_bounds.min[1] = 1;
    global_bounds.min[2] = 0;
    global_bounds.max[0] = 512;
    global_bounds.max[1] = 512;
    global_bounds.max[2] = 0;
    // tells diy to use all availible threads
    //const int num_threads = -1; 
    const int num_threads = 1; 
    const int num_blocks = diy_comm.size(); 
    const int magic_k = 8;

    diy::Master master(diy_comm, num_threads);

    // create an assigner with one block per rank
    diy::ContiguousAssigner assigner(num_blocks, num_blocks); 
    int data = 9595;
    AddIntBlock create(master, data);
    const int num_dims = 2;
    diy::RegularDecomposer<diy::DiscreteBounds> decomposer(num_dims, global_bounds, num_blocks);
    decomposer.decompose(diy_comm.rank(), assigner, create);
    diy::RegularSwapPartners partners(decomposer, 
                                      magic_k, 
                                      false); // false == distance halving
    diy::reduce(master,
                assigner,
                partners,
                debug);

    MPI_Finalize();
}



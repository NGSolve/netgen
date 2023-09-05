#ifndef FILE_ADTREE
#define FILE_ADTREE

/* *************************************************************************/
/* File:   adtree.hh                                                       */
/* Author: Joachim Schoeberl                                               */
/* Date:   16. Feb. 98                                                     */
/* Redesigned by Wolfram Muehlhuber, May 1998                              */
/* *************************************************************************/

#include <general/optmem.hpp>
#include <general/template.hpp>
#include <general/hashtabl.hpp>

#include "geomfuncs.hpp"

namespace netgen
{

/**
  Alternating Digital Tree
 */

// #include "../include/mystdlib.h"
// #include "../include/myadt.hpp"

class ADTreeNode
{
public:
  ADTreeNode *left, *right, *father;
  int dim;
  float sep;
  float *data;
  float *boxmin;
  float *boxmax;
  int pi;
  int nchilds;

  ADTreeNode (int adim);
  ~ADTreeNode ();

  friend class ADTree;
};


class ADTreeCriterion
{
public:
  ADTreeCriterion() { }
  virtual int Eval (const ADTreeNode * node) const = 0;
};


class ADTree
{
  int dim;
  ADTreeNode * root;
  float *cmin, *cmax;
  NgArray<ADTreeNode*> ela;
  const ADTreeCriterion * criterion; 

  NgArray<ADTreeNode*> stack;
  NgArray<int> stackdir;
  int stackindex;

public:
  ADTree (int adim, const float * acmin, 
	   const float * acmax);
  ~ADTree ();

  void Insert (const float * p, int pi);
  // void GetIntersecting (const float * bmin, const float * bmax,
  //			NgArray<int> & pis) const;
  void SetCriterion (ADTreeCriterion & acriterion);
  void Reset ();
  int Next ();
  void GetMatch (NgArray<int> & matches);

  void DeleteElement (int pi);


  void Print (ostream & ost) const
    { PrintRec (ost, root); }

  void PrintRec (ostream & ost, const ADTreeNode * node) const;
};



class ADTreeNode3
{
public:
  ADTreeNode3 *left, *right, *father;
  float sep;
  float data[3];
  int pi;
  int nchilds;

  ADTreeNode3 ();
  void DeleteChilds ();
  friend class ADTree3;

  static BlockAllocator ball;
  void * operator new(size_t);
  void operator delete (void *);
};


class ADTree3
{
  ADTreeNode3 * root;
  float cmin[3], cmax[3];
  NgArray<ADTreeNode3*> ela;

public:
  ADTree3 (const float * acmin, 
	   const float * acmax);
  ~ADTree3 ();

  void Insert (const float * p, int pi);
  void GetIntersecting (const float * bmin, const float * bmax,
			NgArray<int> & pis) const;
  
  void DeleteElement (int pi);


  void Print (ostream & ost) const
    { PrintRec (ost, root); }

  void PrintRec (ostream & ost, const ADTreeNode3 * node) const;
};


/*

// divide each direction
#define ADTN_DIV 10
class ADTreeNode3Div
{
public:
  ADTreeNode3Div *father;
  ADTreeNode3Div *childs[ADTN_DIV];

  float minx, dist;
  float data[3];
  int pi;
  int nchilds;

  ADTreeNode3Div ();
  void DeleteChilds ();
  friend class ADTree3Div;

  static BlockAllocator ball;
  void * operator new(size_t);
  void operator delete (void *);
};


class ADTree3Div
{
  ADTreeNode3Div * root;
  float cmin[3], cmax[3];
  NgArray<ADTreeNode3Div*> ela;

public:
  ADTree3Div (const float * acmin, 
	   const float * acmax);
  ~ADTree3Div ();

  void Insert (const float * p, int pi);
  void GetIntersecting (const float * bmin, const float * bmax,
			NgArray<int> & pis) const;
  
  void DeleteElement (int pi);


  void Print (ostream & ost) const
    { PrintRec (ost, root); }

  void PrintRec (ostream & ost, const ADTreeNode3Div * node) const;
};




#define ADTN_SIZE 10

// multiple entries
class ADTreeNode3M
{
public:
  ADTreeNode3M *left, *right, *father;
  float sep;
  float data[ADTN_SIZE][3];
  int pi[ADTN_SIZE];
  int nchilds;

  ADTreeNode3M ();
  void DeleteChilds ();
  friend class ADTree3M;

  static BlockAllocator ball;
  void * operator new(size_t);
  void operator delete (void *);
};


class ADTree3M
{
  ADTreeNode3M * root;
  float cmin[3], cmax[3];
  NgArray<ADTreeNode3M*> ela;

public:
  ADTree3M (const float * acmin, 
	   const float * acmax);
  ~ADTree3M ();

  void Insert (const float * p, int pi);
  void GetIntersecting (const float * bmin, const float * bmax,
			NgArray<int> & pis) const;
  
  void DeleteElement (int pi);


  void Print (ostream & ost) const
    { PrintRec (ost, root); }

  void PrintRec (ostream & ost, const ADTreeNode3M * node) const;
};






class ADTreeNode3F
{
public:
  ADTreeNode3F *father;
  ADTreeNode3F *childs[8];
  float sep[3];
  float data[3];
  int pi;
  int nchilds;

  ADTreeNode3F ();
  void DeleteChilds ();
  friend class ADTree3F;

  static BlockAllocator ball;
  void * operator new(size_t);
  void operator delete (void *);
};

// fat tree
class ADTree3F
{
  ADTreeNode3F * root;
  float cmin[3], cmax[3];
  NgArray<ADTreeNode3F*> ela;

public:
  ADTree3F (const float * acmin, 
	   const float * acmax);
  ~ADTree3F ();

  void Insert (const float * p, int pi);
  void GetIntersecting (const float * bmin, const float * bmax,
			NgArray<int> & pis) const;
  
  void DeleteElement (int pi);


  void Print (ostream & ost) const
    { PrintRec (ost, root); }

  void PrintRec (ostream & ost, const ADTreeNode3F * node) const;
};




class ADTreeNode3FM
{
public:
  ADTreeNode3FM *father;
  ADTreeNode3FM *childs[8];
  float sep[3];
  float data[ADTN_SIZE][3];
  int pi[ADTN_SIZE];
  int nchilds;

  ADTreeNode3FM ();
  void DeleteChilds ();
  friend class ADTree3FM;

  static BlockAllocator ball;
  void * operator new(size_t);
  void operator delete (void *);
};

// fat tree
class ADTree3FM
{
  ADTreeNode3FM * root;
  float cmin[3], cmax[3];
  NgArray<ADTreeNode3FM*> ela;

public:
  ADTree3FM (const float * acmin, 
	   const float * acmax);
  ~ADTree3FM ();

  void Insert (const float * p, int pi);
  void GetIntersecting (const float * bmin, const float * bmax,
			NgArray<int> & pis) const;
  
  void DeleteElement (int pi);


  void Print (ostream & ost) const
    { PrintRec (ost, root); }

  void PrintRec (ostream & ost, const ADTreeNode3FM * node) const;
};



*/





class ADTreeNode6
{
public:
  ADTreeNode6 *left, *right, *father;
  float sep;
  float data[6];
  int pi;
  int nchilds;

  ADTreeNode6 ();
  void DeleteChilds ();
  friend class ADTree6;

  static BlockAllocator ball;
  void * operator new(size_t);
  void operator delete (void *);
};


class ADTree6
{
  ADTreeNode6 * root;
  float cmin[6], cmax[6];
  NgArray<ADTreeNode6*> ela;

public:
  ADTree6 (const float * acmin, 
	   const float * acmax);
  ~ADTree6 ();

  void Insert (const float * p, int pi);
  void GetIntersecting (const float * bmin, const float * bmax,
			NgArray<int> & pis) const;
  
  void DeleteElement (int pi);

  
  void Print (ostream & ost) const
  { PrintRec (ost, root); }
  int Depth () const
  { return DepthRec (root); }
  int Elements () const
  { return ElementsRec (root); }

  void PrintRec (ostream & ost, const ADTreeNode6 * node) const;
  int DepthRec (const ADTreeNode6 * node) const;
  int ElementsRec (const ADTreeNode6 * node) const;

  void PrintMemInfo (ostream & ost) const;
};





  template <int DIM, typename T>
class T_ADTreeNode
{
public:
  T_ADTreeNode *left, *right, *father;
  float sep;
  // float data[DIM];
  Point<DIM,float> data;
  T pi;
  int nchilds;

  T_ADTreeNode ()
  {
    // pi = -1;
    SetInvalid(pi);
    left = NULL;
    right = NULL;
    father = NULL;
    nchilds = 0;
  }
  
  void DeleteChilds (BlockAllocator & ball)
  {
    if (left)
      {
	left->DeleteChilds(ball);
        ball.Free(left);
	left = NULL;
      }
    if (right)
      {
	right->DeleteChilds(ball);
        ball.Free(right);
	right = NULL;
      }
  }
};



  template <int dim, typename T = INDEX>
  class T_ADTree
  {
    T_ADTreeNode<dim,T> * root;
    // float cmin[dim], cmax[dim];
    Point<dim> cmin, cmax;
    // NgArray<T_ADTreeNode<dim>*> ela;
    NgClosedHashTable<T, T_ADTreeNode<dim,T>*> ela;

    BlockAllocator ball{sizeof(T_ADTreeNode<dim,T>)};
  public:
    T_ADTree (Point<dim> acmin, Point<dim> acmax)
    {
      cmin = acmin;
      cmax = acmax;
      
      root = new (ball.Alloc()) T_ADTreeNode<dim,T>;
      root->sep = (cmin[0] + cmax[0]) / 2;
    }
    
    ~T_ADTree ()
    {
      root->DeleteChilds(ball);
      ball.Free(root);
    }
    
    void Insert (Point<dim> p, T pi)
    {
      T_ADTreeNode<dim,T> *node(NULL);
      T_ADTreeNode<dim,T> *next;
      int dir;
      int lr(0);
      
      Point<dim> bmin = cmin;
      Point<dim> bmax = cmax;
      
      next = root;
      dir = 0;
      while (next)
        {
          node = next;
          
          if (IsInvalid(node->pi))
            {    
              // memcpy (node->data, p, dim * sizeof(float));
              node->data = p;
              node->pi = pi;
              
              // if (ela.Size() < pi+1)
              // ela.SetSize (pi+1);
              ela[pi] = node;
              
	    return;
            }
          
          if (node->sep > p[dir])
            {
              next = node->left;
              bmax(dir) = node->sep;
              lr = 0;
            }
          else
            {
              next = node->right;
              bmin(dir) = node->sep;
              lr = 1;
            }
          
          dir++;
          if (dir == dim) dir = 0;
        }
      
      
      next = new (ball.Alloc()) T_ADTreeNode<dim,T>;
      next->data = p;
      next->pi = pi;
      next->sep = (bmin[dir] + bmax[dir]) / 2;
      
      // if (ela.Size() < pi+1)
      // ela.SetSize (pi+1);
      ela[pi] = next;
      
      if (lr)
        node->right = next;
      else
        node->left = next;
      next -> father = node;

      while (node)
        {
          node->nchilds++;
          node = node->father;
        }
    }

    
    class inttn {
    public:
      int dir;
      T_ADTreeNode<dim,T> * node;
    };
    
    void GetIntersecting (Point<dim> bmin, Point<dim> bmax,
                          NgArray<T> & pis) const
  {
    NgArrayMem<inttn,10000> stack(10000);
    pis.SetSize(0);

    stack[0].node = root;
    stack[0].dir = 0;
    int stacks = 0;

    while (stacks >= 0)
      {
	T_ADTreeNode<dim,T> * node = stack[stacks].node;
	int dir = stack[stacks].dir; 

	stacks--;
	if (!IsInvalid(node->pi)) //  != -1)
	  {
            bool found = true;
            for (int i = 0; i < dim/2; i++)
              if (node->data[i] > bmax[i])
                found = false;
            for (int i = dim/2; i < dim; i++)
              if (node->data[i] < bmin[i])
                found = false;
            if (found)
              pis.Append (node->pi);            
            /*
	    if (node->data[0] > bmax[0] || 
		node->data[1] > bmax[1] || 
		node->data[2] > bmax[2] || 
		node->data[3] < bmin[3] || 
		node->data[4] < bmin[4] || 
		node->data[5] < bmin[5])
	      ;
	    else
              {
                pis.Append (node->pi);
              }
            */
	  }

	int ndir = (dir+1) % dim;

	if (node->left && bmin[dir] <= node->sep)
	  {
	    stacks++;
	    stack[stacks].node = node->left;
	    stack[stacks].dir = ndir;
	  }
	if (node->right && bmax[dir] >= node->sep)
	  {
	    stacks++;
	    stack[stacks].node = node->right;
	    stack[stacks].dir = ndir;
	  }
      }
  }
    
      
    
    void DeleteElement (T pi)
    {
      T_ADTreeNode<dim,T> * node = ela[pi];
      ela.Delete(pi);
      
      SetInvalid(node->pi); //  = -1;
      
      node = node->father;
      while (node)
        {
          node->nchilds--;
          node = node->father;
        }
    }
    
    
    void Print (ostream & ost) const
    { PrintRec (ost, root); }
    int Depth () const
    { return DepthRec (root); }
    int Elements () const
    { return ElementsRec (root); }
    
    void PrintRec (ostream & ost, const T_ADTreeNode<dim,T> * node) const
    {
      
      // if (node->data)     // true anyway
      {
	ost << node->pi << ": ";
	ost << node->nchilds << " childs, ";
	for (int i = 0; i < dim; i++)
	  ost << node->data[i] << " ";
	ost << endl;
      }
      if (node->left)
        PrintRec (ost, node->left);
      if (node->right)
        PrintRec (ost, node->right);
    }

    int DepthRec (const T_ADTreeNode<dim,T> * node) const
    {
      int ldepth = 0;
      int rdepth = 0;
      
      if (node->left)
        ldepth = DepthRec(node->left);
      if (node->right)
        rdepth = DepthRec(node->right);
      return 1 + max2 (ldepth, rdepth);
    }
    
    int ElementsRec (const T_ADTreeNode<dim,T> * node) const
    {
      int els = 1;
      if (node->left)
        els += ElementsRec(node->left);
      if (node->right)
      els += ElementsRec(node->right);
      return els;
    }
    
    
    void PrintMemInfo (ostream & ost) const
    {
      ost << Elements() << " elements a " << sizeof(ADTreeNode6) 
          << " Bytes = "
          << Elements() * sizeof(T_ADTreeNode<dim,T>) << endl;
      ost << "maxind = " << ela.Size() << " = " << sizeof(T_ADTreeNode<dim,T>*) * ela.Size() << " Bytes" << endl;
    }
    
  };
  
  




/*

class ADTreeNode6F
{
public:
  ADTreeNode6F * father;
  ADTreeNode6F * childs[64];
  
  float sep[6];
  float data[6];
  int pi;
  int nchilds;

  ADTreeNode6F ();
  void DeleteChilds ();
  friend class ADTree6F;

  static BlockAllocator ball;
  void * operator new(size_t);
  void operator delete (void *);
};


class ADTree6F
{
  ADTreeNode6F * root;
  float cmin[6], cmax[6];
  NgArray<ADTreeNode6F*> ela;

public:
  ADTree6F (const float * acmin, 
	   const float * acmax);
  ~ADTree6F ();

  void Insert (const float * p, int pi);
  void GetIntersecting (const float * bmin, const float * bmax,
			NgArray<int> & pis) const;
  
  void DeleteElement (int pi);


  void Print (ostream & ost) const
    { PrintRec (ost, root); }
  int Depth () const
    { return DepthRec (root); }

  void PrintRec (ostream & ost, const ADTreeNode6F * node) const;
  int DepthRec (const ADTreeNode6F * node) const;
};







*/





class Point3dTree 
{
  ADTree3 * tree;

public:
  DLL_HEADER Point3dTree (const Point<3> & pmin, const Point<3> & pmax);
  DLL_HEADER ~Point3dTree ();
  DLL_HEADER void Insert (const Point<3> & p, int pi);
  void DeleteElement (int pi) 
    { tree->DeleteElement(pi); }
  DLL_HEADER void GetIntersecting (const Point<3> & pmin, const Point<3> & pmax, 
			NgArray<int> & pis) const;
  const ADTree3 & Tree() const { return *tree; };
};

template<int dim, typename T=INDEX>
class BoxTree
{
public:
  // Number of entries per leaf
  static constexpr int N = 100;

  struct Node;

  struct Leaf
  {
    Point<2*dim> p[N];
    T index[N];
    int n_elements;

    Leaf() : n_elements(0)
    { }

    void Add( NgClosedHashTable<T, Leaf*> &leaf_index, const Point<2*dim> &ap, T aindex )
      {
        p[n_elements] = ap;
        index[n_elements] = aindex;
        n_elements++;
        leaf_index[aindex] = this;
      }
  };

  struct Node
  {
    union
    {
      Node *children[2];
      Leaf *leaf;
    };
    double sep;
    int level;

    Node()
        : children{nullptr,nullptr}
    { }

    ~Node()
     { }

    Leaf *GetLeaf() const
      {
        return children[1] ? nullptr : leaf;
      }
  };

private:
  Node root;

  NgClosedHashTable<T, Leaf*> leaf_index;

  Point<dim> global_min, global_max;
  double tol;
  size_t n_leaves;
  size_t n_nodes;
  BlockAllocator ball_nodes;
  BlockAllocator ball_leaves;

public:

  BoxTree (const Point<dim> & pmin, const Point<dim> & pmax)
      : global_min(pmin), global_max(pmax), n_leaves(1), n_nodes(1), ball_nodes(sizeof(Node)), ball_leaves(sizeof(Leaf))
    {
      root.leaf = (Leaf*) ball_leaves.Alloc(); new (root.leaf) Leaf();
      root.level = 0;
      tol = 1e-7 * Dist(pmax, pmin);
    }

  BoxTree (const Box<dim> & box)
      : BoxTree(box.PMin(), box.PMax())
    { }

  void SetTolerance(double _tol) { tol = _tol; }
  double GetTolerance() { return tol; }

  size_t GetNLeaves()
    {
      return n_leaves;
    }

  size_t GetNNodes()
    {
      return n_nodes;
    }

  template<typename TFunc>
  void GetFirstIntersecting (const Point<dim> & pmin, const Point<dim> & pmax,
          TFunc func=[](auto pi){return false;}) const
    {
      // static Timer timer("BoxTree::GetIntersecting"); RegionTimer rt(timer);
      // static Timer timer1("BoxTree::GetIntersecting-LinearSearch");
      ArrayMem<const Node*, 100> stack;
      ArrayMem<int, 100> dir_stack;


      Point<2*dim> tpmin, tpmax;

      for (size_t i : IntRange(dim))
        {
          tpmin(i) = global_min(i);
          tpmax(i) = pmax(i)+tol;

          tpmin(i+dim) = pmin(i)-tol;
          tpmax(i+dim) = global_max(i);
        }

      stack.SetSize(0);
      stack.Append(&root);
      dir_stack.SetSize(0);
      dir_stack.Append(0);

      while(stack.Size())
        {
          const Node *node = stack.Last();
          stack.DeleteLast();

          int dir = dir_stack.Last();
          dir_stack.DeleteLast();

          if(Leaf *leaf = node->GetLeaf())
            {
              //               RegionTimer rt1(timer1);
              for (auto i : IntRange(leaf->n_elements))
                {
                  bool intersect = true;
                  const auto p = leaf->p[i];

                  for (int d = 0; d < dim; d++)
                      if (p[d] > tpmax[d])
                          intersect = false;
                  for (int d = dim; d < 2*dim; d++)
                      if (p[d] < tpmin[d])
                          intersect = false;
                  if(intersect)
                      if(func(leaf->index[i])) return;
                }
            }
          else
            {
              int newdir = dir+1;
              if(newdir==2*dim) newdir = 0;
              if (tpmin[dir] <= node->sep)
                {
                  stack.Append(node->children[0]);
                  dir_stack.Append(newdir);
                }
              if (tpmax[dir] >= node->sep)
                {
                  stack.Append(node->children[1]);
                  dir_stack.Append(newdir);
                }
            }
        }
    }

  void GetIntersecting (const Point<dim> & pmin, const Point<dim> & pmax,
          NgArray<T> & pis) const
    {
      pis.SetSize(0);
      GetFirstIntersecting(pmin, pmax, [&pis](auto pi) { pis.Append(pi); return false;});
    }

  void GetIntersecting(const Point<dim> & pmin,
                       const Point<dim> & pmax,
                       Array<T> & pis) const
    {
      pis.SetSize0();
      GetFirstIntersecting(pmin, pmax, [&pis](auto pi) { pis.Append(pi); return false;});
    }

  void Insert (const Box<dim> & box, T pi)
    {
      Insert (box.PMin(), box.PMax(), pi);
    }

  void Insert (const Point<dim> & pmin, const Point<dim> & pmax, T pi)
    {
      // static Timer timer("BoxTree::Insert"); RegionTimer rt(timer);
      int dir = 0;
      Point<2*dim> p;
      for (auto i : IntRange(dim))
        {
          p(i) = pmin[i];
          p(i+dim) = pmax[i];
        }

      Node * node = &root;
      Leaf * leaf = node->GetLeaf();

      // search correct leaf to add point
      while(!leaf)
        {
          node = p[dir] < node->sep ? node->children[0] : node->children[1];
          dir++;
          if(dir==2*dim) dir = 0;
          leaf = node->GetLeaf();
        }

      // add point to leaf
      if(leaf->n_elements < N)
          leaf->Add(leaf_index, p,pi);
      else // assume leaf->n_elements == N
        {
          // add two new nodes and one new leaf
          int n_elements = leaf->n_elements;
          ArrayMem<double, N> coords(n_elements);
          ArrayMem<int, N> order(n_elements);

          // separate points in two halves, first sort all coordinates in direction dir
          for (auto i : IntRange(n_elements))
            {
              order[i] = i;
              coords[i] = leaf->p[i][dir];
            }

          QuickSortI(coords, order);
          int isplit = N/2;
          Leaf *leaf1 = (Leaf*) ball_leaves.Alloc(); new (leaf1) Leaf();
          Leaf *leaf2 = (Leaf*) ball_leaves.Alloc(); new (leaf2) Leaf();

          for (auto i : order.Range(isplit))
              leaf1->Add(leaf_index, leaf->p[i], leaf->index[i] );
          for (auto i : order.Range(isplit, N))
              leaf2->Add(leaf_index, leaf->p[i], leaf->index[i] );

          Node *node1 = (Node*) ball_nodes.Alloc(); new (node1) Node();
          node1->leaf = leaf1;
          node1->level = node->level+1;

          Node *node2 = (Node*) ball_nodes.Alloc(); new (node2) Node();
          node2->leaf = leaf2;
          node2->level = node->level+1;

          node->children[0] = node1;
          node->children[1] = node2;
          node->sep = 0.5 * (leaf->p[order[isplit-1]][dir] + leaf->p[order[isplit]][dir]);

          // add new point to one of the new leaves
          if (p[dir] < node->sep)
              leaf1->Add( leaf_index, p, pi );
          else
              leaf2->Add( leaf_index, p, pi );

          ball_leaves.Free(leaf);
          n_leaves++;
          n_nodes+=2;
        }
    }

  void DeleteElement (T pi)
    {
      // static Timer timer("BoxTree::DeleteElement"); RegionTimer rt(timer);
      Leaf *leaf = leaf_index[pi];
      leaf_index.Delete(pi);
      auto & n_elements = leaf->n_elements;
      auto & index = leaf->index;
      auto & p = leaf->p;

      for (auto i : IntRange(n_elements))
        {
          if(index[i] == pi)
            {
              n_elements--;
              if(i!=n_elements)
                {
                  index[i] = index[n_elements];
                  p[i] = p[n_elements];
                }
              return;
            }
        }
    }
};

//   template <int dim, typename T = INDEX>
//   class BoxTree
//   {
//     T_ADTree<2*dim,T> * tree;
//     Point<dim> boxpmin, boxpmax;
//   public:
//     BoxTree (const Box<dim> & abox)
//     {
//       boxpmin = abox.PMin();
//       boxpmax = abox.PMax();
//       Point<2*dim> tpmin, tpmax;
//       for (int i = 0; i < dim; i++)
//         {
//           tpmin(i) = tpmin(i+dim) = boxpmin(i);
//           tpmax(i) = tpmax(i+dim) = boxpmax(i);
//         }
//       tree = new T_ADTree<2*dim,T> (tpmin, tpmax);
//     }
//
//     BoxTree (const Point<dim> & apmin, const Point<dim> & apmax)
//     {
//       boxpmin = apmin;
//       boxpmax = apmax;
//       Point<2*dim> tpmin, tpmax;
//       for (int i = 0; i < dim; i++)
//         {
//           tpmin(i) = tpmin(i+dim) = boxpmin(i);
//           tpmax(i) = tpmax(i+dim) = boxpmax(i);
//         }
//       tree = new T_ADTree<2*dim,T> (tpmin, tpmax);
//     }
//
//     ~BoxTree ()
//     {
//       delete tree;
//     }
//
//     void Insert (const Point<dim> & bmin, const Point<dim> & bmax, T pi)
//     {
//       Point<2*dim> tp;
//
//       for (size_t i = 0; i < dim; i++)
//         {
//           tp(i) = bmin(i);
//           tp(i+dim) = bmax(i);
//         }
//
//       tree->Insert (tp, pi);
//     }
//
//     void Insert (const Box<dim> & box, T pi)
//     {
//       Insert (box.PMin(), box.PMax(), pi);
//     }
//
//     void DeleteElement (T pi)
//     {
//       tree->DeleteElement(pi);
//     }
//
//     void GetIntersecting (const Point<dim> & pmin, const Point<dim> & pmax,
//                           NgArray<T> & pis) const
//     {
//       Point<2*dim> tpmin, tpmax;
//       double tol = Tolerance();
//       for (size_t i = 0; i < dim; i++)
//         {
//           tpmin(i) = boxpmin(i);
//           tpmax(i) = pmax(i)+tol;
//
//           tpmin(i+dim) = pmin(i)-tol;
//           tpmax(i+dim) = boxpmax(i);
//         }
//
//       tree->GetIntersecting (tpmin, tpmax, pis);
//     }
//
//
//     double Tolerance() const { return 1e-7 * Dist(boxpmax, boxpmin); } // single precision
//     const auto & Tree() const { return *tree; };
//     auto & Tree() { return *tree; };
//   };

  template<int dim, typename T=INDEX, typename TSCAL=double>
  class DelaunayTree
  {
  public:
    // Number of entries per leaf
    static constexpr int N = 100;

    struct Node;

    struct Leaf
    {
      Point<2*dim, TSCAL> p[N];
      T index[N];
      int n_elements;
      int nr;

      Leaf() : n_elements(0)
      { }


      void Add( Array<Leaf*> &leaves, Array<T> &leaf_index, const Point<2*dim> &ap, T aindex )
        {
          p[n_elements] = ap;
          index[n_elements] = aindex;
          n_elements++;
          if(leaf_index.Size()<aindex+1)
              leaf_index.SetSize(aindex+1);
          leaf_index[aindex] = nr;
        }
    };

    struct Node
    {
      union
      {
        Node *children[2];
        Leaf *leaf;
      };
      double sep;
      int level;

      Node()
          : children{nullptr,nullptr}
      { }

      ~Node()
       { }

      Leaf *GetLeaf() const
        {
          return children[1] ? nullptr : leaf;
        }
    };

  private:
    Node root;

    Array<Leaf*> leaves;
    Array<T> leaf_index;

    Point<dim> global_min, global_max;
    double tol;
    size_t n_leaves;
    size_t n_nodes;
    BlockAllocator ball_nodes;
    BlockAllocator ball_leaves;

  public:

    DelaunayTree (const Point<dim> & pmin, const Point<dim> & pmax)
        : global_min(pmin), global_max(pmax), n_leaves(1), n_nodes(1), ball_nodes(sizeof(Node)), ball_leaves(sizeof(Leaf))
      {
        root.leaf = (Leaf*) ball_leaves.Alloc(); new (root.leaf) Leaf();
        root.leaf->nr = 0;
        leaves.Append(root.leaf);
        root.level = 0;
        tol = 1e-7 * Dist(pmax, pmin);
      }

    DelaunayTree (const Box<dim> & box)
        : DelaunayTree(box.PMin(), box.PMax())
      { }

    double GetTolerance() { return tol; }

    size_t GetNLeaves()
      {
        return n_leaves;
      }

    size_t GetNNodes()
      {
        return n_nodes;
      }

    template<typename TFunc>
    void GetFirstIntersecting (const Point<dim> & pmin, const Point<dim> & pmax,
            TFunc func=[](auto pi){return false;}) const
      {
        // static Timer timer("DelaunayTree::GetIntersecting"); RegionTimer rt(timer);
        // static Timer timer1("DelaunayTree::GetIntersecting-LinearSearch");
        ArrayMem<const Node*, 100> stack;
        ArrayMem<int, 100> dir_stack;


        Point<2*dim> tpmin, tpmax;

        for (size_t i : IntRange(dim))
          {
            tpmin(i) = global_min(i);
            tpmax(i) = pmax(i)+tol;

            tpmin(i+dim) = pmin(i)-tol;
            tpmax(i+dim) = global_max(i);
          }

        stack.SetSize(0);
        stack.Append(&root);
        dir_stack.SetSize(0);
        dir_stack.Append(0);

        while(stack.Size())
          {
            const Node *node = stack.Last();
            stack.DeleteLast();

            int dir = dir_stack.Last();
            dir_stack.DeleteLast();

            if(Leaf *leaf = node->GetLeaf())
              {
                //               RegionTimer rt1(timer1);
                for (auto i : IntRange(leaf->n_elements))
                  {
                    bool intersect = true;
                    const auto p = leaf->p[i];

                    for (int d = 0; d < dim; d++)
                        if (p[d] > tpmax[d])
                            intersect = false;
                    for (int d = dim; d < 2*dim; d++)
                        if (p[d] < tpmin[d])
                            intersect = false;
                    if(intersect)
                        if(func(leaf->index[i])) return;
                  }
              }
            else
              {
                int newdir = dir+1;
                if(newdir==2*dim) newdir = 0;
                if (tpmin[dir] <= node->sep)
                  {
                    stack.Append(node->children[0]);
                    dir_stack.Append(newdir);
                  }
                if (tpmax[dir] >= node->sep)
                  {
                    stack.Append(node->children[1]);
                    dir_stack.Append(newdir);
                  }
              }
          }
      }

    void GetIntersecting (const Point<dim> & pmin, const Point<dim> & pmax,
            NgArray<T> & pis) const
      {
        pis.SetSize(0);
        GetFirstIntersecting(pmin, pmax, [&pis](auto pi) { pis.Append(pi); return false;});
      }

    void Insert (const Box<dim> & box, T pi)
      {
        Insert (box.PMin(), box.PMax(), pi);
      }

    void Insert (const Point<dim> & pmin, const Point<dim> & pmax, T pi)
      {
        // static Timer timer("DelaunayTree::Insert"); RegionTimer rt(timer);
        int dir = 0;
        Point<2*dim> p;
        for (auto i : IntRange(dim))
          {
            p(i) = pmin[i];
            p(i+dim) = pmax[i];
          }

        Node * node = &root;
        Leaf * leaf = node->GetLeaf();

        // search correct leaf to add point
        while(!leaf)
          {
            node = p[dir] < node->sep ? node->children[0] : node->children[1];
            dir++;
            if(dir==2*dim) dir = 0;
            leaf = node->GetLeaf();
          }

        // add point to leaf
        if(leaf->n_elements < N)
            leaf->Add(leaves, leaf_index, p,pi);
        else // assume leaf->n_elements == N
          {
            // add two new nodes and one new leaf
            int n_elements = leaf->n_elements;
            ArrayMem<TSCAL, N> coords(n_elements);
            ArrayMem<int, N> order(n_elements);

            // separate points in two halves, first sort all coordinates in direction dir
            for (auto i : IntRange(n_elements))
              {
                order[i] = i;
                coords[i] = leaf->p[i][dir];
              }

            QuickSortI(coords, order);
            int isplit = N/2;
            Leaf *leaf1 = (Leaf*) ball_leaves.Alloc(); new (leaf1) Leaf();
            Leaf *leaf2 = (Leaf*) ball_leaves.Alloc(); new (leaf2) Leaf();

            leaf1->nr = leaf->nr;
            leaf2->nr = leaves.Size();
            leaves.Append(leaf2);
            leaves[leaf1->nr] = leaf1;

            for (auto i : order.Range(isplit))
                leaf1->Add(leaves, leaf_index, leaf->p[i], leaf->index[i] );
            for (auto i : order.Range(isplit, N))
                leaf2->Add(leaves, leaf_index, leaf->p[i], leaf->index[i] );

            Node *node1 = (Node*) ball_nodes.Alloc(); new (node1) Node();
            node1->leaf = leaf1;
            node1->level = node->level+1;

            Node *node2 = (Node*) ball_nodes.Alloc(); new (node2) Node();
            node2->leaf = leaf2;
            node2->level = node->level+1;

            node->children[0] = node1;
            node->children[1] = node2;
            node->sep = 0.5 * (leaf->p[order[isplit-1]][dir] + leaf->p[order[isplit]][dir]);

            // add new point to one of the new leaves
            if (p[dir] < node->sep)
                leaf1->Add( leaves, leaf_index, p, pi );
            else
                leaf2->Add( leaves, leaf_index, p, pi );

            ball_leaves.Free(leaf);
            n_leaves++;
            n_nodes+=2;
          }
      }

    void DeleteElement (T pi)
      {
        // static Timer timer("DelaunayTree::DeleteElement"); RegionTimer rt(timer);
        Leaf *leaf = leaves[leaf_index[pi]];
        leaf_index[pi] = -1;
        auto & n_elements = leaf->n_elements;
        auto & index = leaf->index;
        auto & p = leaf->p;

        for (auto i : IntRange(n_elements))
          {
            if(index[i] == pi)
              {
                n_elements--;
                if(i!=n_elements)
                  {
                    index[i] = index[n_elements];
                    p[i] = p[n_elements];
                  }
                return;
              }
          }
      }
  };
}

#endif

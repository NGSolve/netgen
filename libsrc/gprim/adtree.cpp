#include <mystdlib.h>


#include <myadt.hpp>
// class DenseMatrix;
#include <gprim.hpp>

namespace netgen
{


  /* ******************************* ADTree ******************************* */


  ADTreeNode :: ADTreeNode(int adim)
  {
    pi = -1;

    left = NULL;
    right = NULL;
    father = NULL;
    nchilds = 0;
    dim = adim;
    data = new float [dim];
    boxmin = NULL;
    boxmax = NULL;
  }




  ADTreeNode :: ~ADTreeNode()
  {
    delete data;
  }


  ADTree :: ADTree (int adim, const float * acmin, 
		    const float * acmax)
    : ela(0), stack(1000), stackdir(1000)
  {
    dim = adim;
    cmin = new float [dim];
    cmax = new float [dim];
    memcpy (cmin, acmin, dim * sizeof(float));
    memcpy (cmax, acmax, dim * sizeof(float));

    root = new ADTreeNode (dim);
    root->sep = (cmin[0] + cmax[0]) / 2;
    root->boxmin = new float [dim];
    root->boxmax = new float [dim];
    memcpy (root->boxmin, cmin, dim * sizeof(float));
    memcpy (root->boxmax, cmax, dim * sizeof(float));
  }

  ADTree :: ~ADTree ()
  {
    ;
  }

  void ADTree :: Insert (const float * p, int pi)
  {
    ADTreeNode *node(NULL);
    ADTreeNode *next;
    int dir;
    int lr(1);

    float * bmin = new float [dim];
    float * bmax = new float [dim];
  
    memcpy (bmin, cmin, dim * sizeof(float));
    memcpy (bmax, cmax, dim * sizeof(float));


    next = root;
    dir = 0;
    while (next)
      {
	node = next;

	if (node->pi == -1)
	  {    
	    memcpy (node->data, p, dim * sizeof(float));
	    node->pi = pi;

	    if (ela.Size() < pi+1)
	      ela.SetSize (pi+1);
	    ela[pi] = node;

	    return;
	  }

	if (node->sep > p[dir])
	  {
	    next = node->left;
	    bmax[dir] = node->sep;
	    lr = 0;
	  }
	else
	  {
	    next = node->right;
	    bmin[dir] = node->sep;
	    lr = 1;
	  }

	dir++;
	if (dir == dim)
	  dir = 0;
      }


    next = new ADTreeNode(dim);
    memcpy (next->data, p, dim * sizeof(float));
    next->pi = pi;
    next->sep = (bmin[dir] + bmax[dir]) / 2;
    next->boxmin = bmin;
    next->boxmax = bmax;

    if (ela.Size() < pi+1)
      ela.SetSize (pi+1);
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

  void ADTree :: DeleteElement (int pi)
  {
    ADTreeNode * node = ela[pi];

    node->pi = -1;

    node = node->father;
    while (node)
      {
	node->nchilds--;
	node = node->father;
      }
  }


  void ADTree :: SetCriterion (ADTreeCriterion & acriterion)
  {
    criterion = & acriterion;
  }


  void ADTree :: Reset ()
  {
    stack.Elem(1) = root;
    stackdir.Elem(1) = 0;
    stackindex = 1;
  }


  int ADTree:: Next ()
  {
    ADTreeNode *node;
    int dir;

    if (stackindex == 0)
      return -1;

    do 
      {
	node = stack.Get(stackindex);
	dir = stackdir.Get(stackindex);
	stackindex --;

	if (criterion -> Eval(node))
	  {
	    int ndir = dir + 1;
	    if (ndir == dim)
	      ndir = 0;

	    if (node -> left && criterion -> Eval (node->left))
	      {
		stackindex ++;
		stack.Elem(stackindex) = node -> left;
		stackdir.Elem(stackindex) = ndir;
	      }
	    if (node->right && criterion -> Eval (node -> right))
	      {
		stackindex++;
		stack.Elem(stackindex) = node->right;
		stackdir.Elem(stackindex) = ndir;
	      }
	  
	    if (node -> pi != -1)
	      return node->pi;
	  }
      }
    while (stackindex > 0);

    return -1;
  }


  void ADTree :: GetMatch (NgArray <int> & matches)
  {
    int nodenr;

    Reset();

    while ( (nodenr = Next()) != -1)
      matches.Append (nodenr);
  }


  void ADTree :: PrintRec (ostream & ost, const ADTreeNode * node) const
  {
  
    if (node->data)
      {
	ost << node->pi << ": ";
	ost << node->nchilds << " childs, ";
	for (int i = 0; i < dim; i++)
	  ost << node->data[i] << " ";
	ost << endl;
      }
    if (node->left)
      {
	ost << "l ";
	PrintRec (ost, node->left);
      }
    if (node->right)
      {
	ost << "r ";
	PrintRec (ost, node->right);
      }
  }


  /* ******************************* ADTree3 ******************************* */


  ADTreeNode3 :: ADTreeNode3()
  {
    pi = -1;

    left = NULL;
    right = NULL;
    father = NULL;
    nchilds = 0;
  }

  void ADTreeNode3 :: DeleteChilds ()
  {
    if (left)
      {
	left->DeleteChilds();
	delete left;
	left = NULL;
      }
    if (right)
      {
	right->DeleteChilds();
	delete right;
	right = NULL;
      }
  }


  BlockAllocator ADTreeNode3 :: ball(sizeof (ADTreeNode3));


  void * ADTreeNode3 :: operator new(size_t s)
  {
    return ball.Alloc();
  }

  void ADTreeNode3 :: operator delete (void * p)
  {
    ball.Free (p);
  }







  ADTree3 :: ADTree3 (const float * acmin, 
		      const float * acmax)
    : ela(0)
  {
    memcpy (cmin, acmin, 3 * sizeof(float));
    memcpy (cmax, acmax, 3 * sizeof(float));

    root = new ADTreeNode3;
    root->sep = (cmin[0] + cmax[0]) / 2;
  }

  ADTree3 :: ~ADTree3 ()
  {
    root->DeleteChilds();
    delete root;
  }


  void ADTree3 :: Insert (const float * p, int pi)
  {
    ADTreeNode3 *node(NULL);
    ADTreeNode3 *next;
    int dir;
    int lr(0);

    float bmin[3];
    float bmax[3];
  
    memcpy (bmin, cmin, 3 * sizeof(float));
    memcpy (bmax, cmax, 3 * sizeof(float));

    next = root;
    dir = 0;
    while (next)
      {
	node = next;

	if (node->pi == -1)
	  {    
	    memcpy (node->data, p, 3 * sizeof(float));
	    node->pi = pi;

	    if (ela.Size() < pi+1)
	      ela.SetSize (pi+1);
	    ela[pi] = node;

	    return;
	  }

	if (node->sep > p[dir])
	  {
	    next = node->left;
	    bmax[dir] = node->sep;
	    lr = 0;
	  }
	else
	  {
	    next = node->right;
	    bmin[dir] = node->sep;
	    lr = 1;
	  }

	dir++;
	if (dir == 3)
	  dir = 0;
      }


    next = new ADTreeNode3;
    memcpy (next->data, p, 3 * sizeof(float));
    next->pi = pi;
    next->sep = (bmin[dir] + bmax[dir]) / 2;


    if (ela.Size() < pi+1)
      ela.SetSize (pi+1);
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

  void ADTree3 :: DeleteElement (int pi)
  {
    ADTreeNode3 * node = ela[pi];

    node->pi = -1;

    node = node->father;
    while (node)
      {
	node->nchilds--;
	node = node->father;
      }
  }

  void ADTree3 :: GetIntersecting (const float * bmin, 
				   const float * bmax,
				   NgArray<int> & pis) const
  {
    ArrayMem<ADTreeNode3*, 1000> stack(1000);
    ArrayMem<int, 1000> stackdir(1000);
    ADTreeNode3 * node;
    int dir, stacks;

    // stack.SetSize (1000);
    // stackdir.SetSize(1000);
    pis.SetSize(0);

    stack[0] = root;
    stackdir[0] = 0;
    stacks = 0;

    while (stacks >= 0)
      {
	node = stack[stacks];
	dir = stackdir[stacks]; 
	stacks--;

	if (node->pi != -1)
	  {
	    if (node->data[0] >= bmin[0] && node->data[0] <= bmax[0] &&
		node->data[1] >= bmin[1] && node->data[1] <= bmax[1] &&
		node->data[2] >= bmin[2] && node->data[2] <= bmax[2])

	      pis.Append (node->pi);
	  }


	int ndir = dir+1;
	if (ndir == 3)
	  ndir = 0;

	if (node->left && bmin[dir] <= node->sep)
	  {
	    stacks++;
	    stack[stacks] = node->left;
	    stackdir[stacks] = ndir;
	  }
	if (node->right && bmax[dir] >= node->sep)
	  {
	    stacks++;
	    stack[stacks] = node->right;
	    stackdir[stacks] = ndir;
	  }
      }
  }

  void ADTree3 :: PrintRec (ostream & ost, const ADTreeNode3 * node) const
  {
  
    // if (node->data)  // true anyway
      {
	ost << node->pi << ": ";
	ost << node->nchilds << " childs, ";
	for (int i = 0; i < 3; i++)
	  ost << node->data[i] << " ";
	ost << endl;
      }
    if (node->left)
      PrintRec (ost, node->left);
    if (node->right)
      PrintRec (ost, node->right);
  }








#ifdef ABC

  /* ******************************* ADTree3Div ******************************* */


  ADTreeNode3Div :: ADTreeNode3Div()
  {
    pi = 0;
  
    int i;
    for (i = 0; i < ADTN_DIV; i++)
      childs[i] = NULL;

    father = NULL;
    nchilds = 0;
    minx = 0;
    dist = 1;
  }

  void ADTreeNode3Div :: DeleteChilds ()
  {
    int i;
    for (i = 0; i < ADTN_DIV; i++)
      if (childs[i])
	{
	  childs[i]->DeleteChilds();
	  delete childs[i];
	  childs[i] = NULL;
	}
  }


  BlockAllocator ADTreeNode3Div :: ball(sizeof (ADTreeNode3Div));

  void * ADTreeNode3Div :: operator new(size_t)
  {
    return ball.Alloc();
  }

  void ADTreeNode3Div :: operator delete (void * p)
  {
    ball.Free (p);
  }







  ADTree3Div :: ADTree3Div (const float * acmin, 
			    const float * acmax)
    : ela(0)
  {
    memcpy (cmin, acmin, 3 * sizeof(float));
    memcpy (cmax, acmax, 3 * sizeof(float));

    root = new ADTreeNode3Div;

    root->minx = cmin[0];
    root->dist = (cmax[0] - cmin[0]) / ADTN_DIV;

    //  root->sep = (cmin[0] + cmax[0]) / 2;
  }

  ADTree3Div :: ~ADTree3Div ()
  {
    root->DeleteChilds();
    delete root;
  }


  void ADTree3Div :: Insert (const float * p, int pi)
  {
    ADTreeNode3Div *node;
    ADTreeNode3Div *next;
    int dir;
    int bag;
  
    float bmin[3];
    float bmax[3];
  
    memcpy (bmin, cmin, 3 * sizeof(float));
    memcpy (bmax, cmax, 3 * sizeof(float));


    next = root;
    dir = 0;
    while (next)
      {
	node = next;

	if (!node->pi)
	  {    
	    memcpy (node->data, p, 3 * sizeof(float));
	    node->pi = pi;

	    if (ela.Size() < pi)
	      ela.SetSize (pi);
	    ela.Elem(pi) = node;

	    return;
	  }

	double dx = (bmax[dir] - bmin[dir]) / ADTN_DIV;
	bag = int ((p[dir]-bmin[dir]) / dx);

	//      (*testout) << "insert, bag = " << bag << endl;

	if (bag < 0) bag = 0;
	if (bag >= ADTN_DIV) bag = ADTN_DIV-1;
      
	double nbmin = bmin[dir] + bag * dx;
	double nbmax = bmin[dir] + (bag+1) * dx;

	/*      
		(*testout) << "bmin, max = " << bmin[dir] << "-" << bmax[dir]
		<< " p = " << p[dir];
	*/
	next = node->childs[bag];
	bmin[dir] = nbmin;
	bmax[dir] = nbmax;

	//      (*testout) << "new bmin, max = " << bmin[dir] << "-" << bmax[dir] << endl;

      
	/*      
		if (node->sep > p[dir])
		{
		next = node->left;
		bmax[dir] = node->sep;
		lr = 0;
		}
		else
		{
		next = node->right;
		bmin[dir] = node->sep;
		lr = 1;
		}
	*/

	dir++;
	if (dir == 3)
	  dir = 0;
      }


    next = new ADTreeNode3Div;
    memcpy (next->data, p, 3 * sizeof(float));
    next->pi = pi;

    next->minx = bmin[dir];
    next->dist = (bmax[dir] - bmin[dir]) / ADTN_DIV;
    //  next->sep = (bmin[dir] + bmax[dir]) / 2;


    if (ela.Size() < pi)
      ela.SetSize (pi);
    ela.Elem(pi) = next;

    node->childs[bag] = next;
    next -> father = node;

    while (node)
      {
	node->nchilds++;
	node = node->father;
      }
  }

  void ADTree3Div :: DeleteElement (int pi)
  {
    ADTreeNode3Div * node = ela.Get(pi);

    node->pi = 0;

    node = node->father;
    while (node)
      {
	node->nchilds--;
	node = node->father;
      }
  }

  void ADTree3Div :: GetIntersecting (const float * bmin, 
				      const float * bmax,
				      NgArray<int> & pis) const
  {
    static NgArray<ADTreeNode3Div*> stack(1000);
    static NgArray<int> stackdir(1000);
    ADTreeNode3Div * node;
    int dir, i, stacks;

    stack.SetSize (1000);
    stackdir.SetSize(1000);
    pis.SetSize(0);

    stack.Elem(1) = root;
    stackdir.Elem(1) = 0;
    stacks = 1;

    while (stacks)
      {
	node = stack.Get(stacks);
	dir = stackdir.Get(stacks); 
	stacks--;

	if (node->pi)
	  {
	    if (node->data[0] >= bmin[0] && node->data[0] <= bmax[0] &&
		node->data[1] >= bmin[1] && node->data[1] <= bmax[1] &&
		node->data[2] >= bmin[2] && node->data[2] <= bmax[2])

	      pis.Append (node->pi);
	  }


	int ndir = dir+1;
	if (ndir == 3)
	  ndir = 0;

	int mini = int ( (bmin[dir] - node->minx) / node->dist );
	int maxi = int ( (bmax[dir] - node->minx) / node->dist );
      
	//      (*testout) << "get int, mini, maxi = " << mini << ", " << maxi << endl;
	if (mini < 0) mini = 0;
	if (maxi >= ADTN_DIV) maxi = ADTN_DIV-1;

	for (i = mini; i <= maxi; i++)
	  if (node->childs[i])
	    {
	      stacks++;
	      stack.Elem(stacks) = node->childs[i];
	      stackdir.Elem(stacks) = ndir;
	    }


	/*
	  if (node->left && bmin[dir] <= node->sep)
	  {
	  stacks++;
	  stack.Elem(stacks) = node->left;
	  stackdir.Elem(stacks) = ndir;
	  }
	  if (node->right && bmax[dir] >= node->sep)
	  {
	  stacks++;
	  stack.Elem(stacks) = node->right;
	  stackdir.Elem(stacks) = ndir;
	  }
	*/
      }
  }

  void ADTree3Div :: PrintRec (ostream & ost, const ADTreeNode3Div * node) const
  {
  
    if (node->data)
      {
	ost << node->pi << ": ";
	ost << node->nchilds << " childs, ";
	ost << " from " << node->minx << " - " << node->minx + node->dist*ADTN_DIV << "  ";
	for (int i = 0; i < 3; i++)
	  ost << node->data[i] << " ";
	ost << endl;
      }
    int i;
    for (i = 0; i < ADTN_DIV; i++)
      if (node->childs[i])
	PrintRec (ost, node->childs[i]);
  }












  /* ******************************* ADTree3M ******************************* */


  ADTreeNode3M :: ADTreeNode3M()
  {
    int i;
    for (i = 0; i < ADTN_SIZE; i++)
      pi[i] = 0;

    left = NULL;
    right = NULL;
    father = NULL;
    nchilds = 0;
  }

  void ADTreeNode3M :: DeleteChilds ()
  {
    if (left)
      {
	left->DeleteChilds();
	delete left;
	left = NULL;
      }
    if (right)
      {
	right->DeleteChilds();
	delete right;
	right = NULL;
      }
  }


  BlockAllocator ADTreeNode3M :: ball(sizeof (ADTreeNode3M));

  void * ADTreeNode3M :: operator new(size_t)
  {
    return ball.Alloc();
  }

  void ADTreeNode3M :: operator delete (void * p)
  {
    ball.Free (p);
  }







  ADTree3M :: ADTree3M (const float * acmin, 
			const float * acmax)
    : ela(0)
  {
    memcpy (cmin, acmin, 3 * sizeof(float));
    memcpy (cmax, acmax, 3 * sizeof(float));

    root = new ADTreeNode3M;
    root->sep = (cmin[0] + cmax[0]) / 2;
  }

  ADTree3M :: ~ADTree3M ()
  {
    root->DeleteChilds();
    delete root;
  }


  void ADTree3M :: Insert (const float * p, int pi)
  {
    ADTreeNode3M *node;
    ADTreeNode3M *next;
    int dir;
    int lr;
    int i;
    float bmin[3];
    float bmax[3];
  
    memcpy (bmin, cmin, 3 * sizeof(float));
    memcpy (bmax, cmax, 3 * sizeof(float));

    next = root;
    dir = 0;
    while (next)
      {
	node = next;

	for (i = 0; i < ADTN_SIZE; i++)
	  if (!node->pi[i])
	    {    
	      memcpy (node->data[i], p, 3 * sizeof(float));
	      node->pi[i] = pi;
	    
	      if (ela.Size() < pi)
		ela.SetSize (pi);
	      ela.Elem(pi) = node;
	    
	      return;
	    }

	if (node->sep > p[dir])
	  {
	    next = node->left;
	    bmax[dir] = node->sep;
	    lr = 0;
	  }
	else
	  {
	    next = node->right;
	    bmin[dir] = node->sep;
	    lr = 1;
	  }

	dir++;
	if (dir == 3)
	  dir = 0;
      }


    next = new ADTreeNode3M;
    memcpy (next->data[0], p, 3 * sizeof(float));
    next->pi[0] = pi;
    next->sep = (bmin[dir] + bmax[dir]) / 2;


    if (ela.Size() < pi)
      ela.SetSize (pi);
    ela.Elem(pi) = next;


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

  void ADTree3M :: DeleteElement (int pi)
  {
    ADTreeNode3M * node = ela.Get(pi);

    int i;
    for (i = 0; i < ADTN_SIZE; i++)
      if (node->pi[i] == pi)
	node->pi[i] = 0;

    node = node->father;
    while (node)
      {
	node->nchilds--;
	node = node->father;
      }
  }

  void ADTree3M :: GetIntersecting (const float * bmin, 
				    const float * bmax,
				    NgArray<int> & pis) const
  {
    static NgArray<ADTreeNode3M*> stack(1000);
    static NgArray<int> stackdir(1000);
    ADTreeNode3M * node;
    int dir, i, stacks;

    stack.SetSize (1000);
    stackdir.SetSize(1000);
    pis.SetSize(0);

    stack.Elem(1) = root;
    stackdir.Elem(1) = 0;
    stacks = 1;

    while (stacks)
      {
	node = stack.Get(stacks);
	dir = stackdir.Get(stacks); 
	stacks--;

	int * hpi = node->pi;
	for (i = 0; i < ADTN_SIZE; i++)
	  if (hpi[i])
	    {
	      float * datai = &node->data[i][0];
	      if (datai[0] >= bmin[0] && datai[0] <= bmax[0] &&
		  datai[1] >= bmin[1] && datai[1] <= bmax[1] &&
		  datai[2] >= bmin[2] && datai[2] <= bmax[2])
	      
		pis.Append (node->pi[i]);
	    }


	int ndir = dir+1;
	if (ndir == 3)
	  ndir = 0;

	if (node->left && bmin[dir] <= node->sep)
	  {
	    stacks++;
	    stack.Elem(stacks) = node->left;
	    stackdir.Elem(stacks) = ndir;
	  }
	if (node->right && bmax[dir] >= node->sep)
	  {
	    stacks++;
	    stack.Elem(stacks) = node->right;
	    stackdir.Elem(stacks) = ndir;
	  }
      }
  }

  void ADTree3M :: PrintRec (ostream & ost, const ADTreeNode3M * node) const
  {
  
    if (node->data)
      {
	//      ost << node->pi << ": ";
	ost << node->nchilds << " childs, ";
	for (int i = 0; i < 3; i++)
	  ost << node->data[i] << " ";
	ost << endl;
      }
    if (node->left)
      PrintRec (ost, node->left);
    if (node->right)
      PrintRec (ost, node->right);
  }












  /* ******************************* ADTree3F ******************************* */


  ADTreeNode3F :: ADTreeNode3F()
  {
    pi = 0;
    father = NULL;
    nchilds = 0;
    int i;
    for (i = 0; i < 8; i++)
      childs[i] = NULL;
  }

  void ADTreeNode3F :: DeleteChilds ()
  {
    int i;

    for (i = 0; i < 8; i++)
      {
	if (childs[i])
	  childs[i]->DeleteChilds();
	delete childs[i];
	childs[i] = NULL;
      }
  }


  BlockAllocator ADTreeNode3F :: ball(sizeof (ADTreeNode3F));

  void * ADTreeNode3F :: operator new(size_t)
  {
    return ball.Alloc();
  }

  void ADTreeNode3F :: operator delete (void * p)
  {
    ball.Free (p);
  }







  ADTree3F :: ADTree3F (const float * acmin, 
			const float * acmax)
    : ela(0)
  {
    memcpy (cmin, acmin, 3 * sizeof(float));
    memcpy (cmax, acmax, 3 * sizeof(float));

    root = new ADTreeNode3F;
    for (int i = 0; i < 3; i++)
      root->sep[i] = (cmin[i] + cmax[i]) / 2;
  }

  ADTree3F :: ~ADTree3F ()
  {
    root->DeleteChilds();
    delete root;
  }


  void ADTree3F :: Insert (const float * p, int pi)
  {
    ADTreeNode3F *node;
    ADTreeNode3F *next;
    int lr;

    float bmin[3];
    float bmax[3];
    int i, dir;
  
    memcpy (bmin, cmin, 3 * sizeof(float));
    memcpy (bmax, cmax, 3 * sizeof(float));


    next = root;
    while (next)
      {
	node = next;
      
	if (!node->pi)
	  {    
	    memcpy (node->data, p, 3 * sizeof(float));
	    node->pi = pi;

	    if (ela.Size() < pi)
	      ela.SetSize (pi);
	    ela.Elem(pi) = node;

	    return;
	  }

	dir = 0;
	for (i = 0; i < 3; i++)
	  {
	    if (node->sep[i] > p[i])
	      {
		bmax[i] = node->sep[i];
	      }
	    else
	      {
		bmin[i] = node->sep[i];
		dir += (1 << i);
	      }
	  }
	next = node->childs[dir];

	/*
	  if (node->sep > p[dir])
	  {
	  next = node->left;
	  bmax[dir] = node->sep;
	  lr = 0;
	  }
	  else
	  {
	  next = node->right;
	  bmin[dir] = node->sep;
	  lr = 1;
	  }
	*/
      }


    next = new ADTreeNode3F;
    memcpy (next->data, p, 3 * sizeof(float));
    next->pi = pi;

    for (i = 0; i < 3; i++)
      next->sep[i] = (bmin[i] + bmax[i]) / 2;
  

    if (ela.Size() < pi)
      ela.SetSize (pi);
    ela.Elem(pi) = next;

    node->childs[dir] = next;
    next->father = node;

    while (node)
      {
	node->nchilds++;
	node = node->father;
      }
  }

  void ADTree3F :: DeleteElement (int pi)
  {
    ADTreeNode3F * node = ela.Get(pi);

    node->pi = 0;

    node = node->father;
    while (node)
      {
	node->nchilds--;
	node = node->father;
      }
  }

  void ADTree3F :: GetIntersecting (const float * bmin, 
				    const float * bmax,
				    NgArray<int> & pis) const
  {
    static NgArray<ADTreeNode3F*> stack(1000);
    ADTreeNode3F * node;
    int dir, i, stacks;

    stack.SetSize (1000);
    pis.SetSize(0);

    stack.Elem(1) = root;
    stacks = 1;

    while (stacks)
      {
	node = stack.Get(stacks);
	stacks--;

	if (node->pi)
	  {
	    if (node->data[0] >= bmin[0] && node->data[0] <= bmax[0] &&
		node->data[1] >= bmin[1] && node->data[1] <= bmax[1] &&
		node->data[2] >= bmin[2] && node->data[2] <= bmax[2])

	      pis.Append (node->pi);
	  }

      
	int i1min = (bmin[0] <= node->sep[0]) ? 0 : 1;
	int i1max = (bmax[0] < node->sep[0]) ? 0 : 1;
	int i2min = (bmin[1] <= node->sep[1]) ? 0 : 1;
	int i2max = (bmax[1] < node->sep[1]) ? 0 : 1;
	int i3min = (bmin[2] <= node->sep[2]) ? 0 : 1;
	int i3max = (bmax[2] < node->sep[2]) ? 0 : 1;

	int i1, i2, i3;
	for (i1 = i1min; i1 <= i1max; i1++)
	  for (i2 = i2min; i2 <= i2max; i2++)
	    for (i3 = i3min; i3 <= i3max; i3++)
	      {
		i = i1+2*i2+4*i3;
		if (node->childs[i])
		  {
		    stacks++;
		    stack.Elem(stacks) = node->childs[i];
		  }
	      }
      
	/*
	  if (node->left && bmin[dir] <= node->sep)
	  {
	  stacks++;
	  stack.Elem(stacks) = node->left;
	  stackdir.Elem(stacks) = ndir;
	  }
	  if (node->right && bmax[dir] >= node->sep)
	  {
	  stacks++;
	  stack.Elem(stacks) = node->right;
	  stackdir.Elem(stacks) = ndir;
	  }
	*/
      }
  }

  void ADTree3F :: PrintRec (ostream & ost, const ADTreeNode3F * node) const
  {
    int i;
    if (node->data)
      {
	ost << node->pi << ": ";
	ost << node->nchilds << " childs, ";
	for (i = 0; i < 3; i++)
	  ost << node->data[i] << " ";
	ost << endl;
      }

    for (i = 0; i < 8; i++)
      if (node->childs[i])
	PrintRec (ost, node->childs[i]);
  }













  /* ******************************* ADTree3FM ******************************* */


  ADTreeNode3FM :: ADTreeNode3FM()
  {
    father = NULL;
    nchilds = 0;
    int i;

    for (i = 0; i < ADTN_SIZE; i++)
      pi[i] = 0;

    for (i = 0; i < 8; i++)
      childs[i] = NULL;
  }

  void ADTreeNode3FM :: DeleteChilds ()
  {
    int i;

    for (i = 0; i < 8; i++)
      {
	if (childs[i])
	  childs[i]->DeleteChilds();
	delete childs[i];
	childs[i] = NULL;
      }
  }


  BlockAllocator ADTreeNode3FM :: ball(sizeof (ADTreeNode3FM));

  void * ADTreeNode3FM :: operator new(size_t)
  {
    return ball.Alloc();
  }

  void ADTreeNode3FM :: operator delete (void * p)
  {
    ball.Free (p);
  }







  ADTree3FM :: ADTree3FM (const float * acmin, 
			  const float * acmax)
    : ela(0)
  {
    memcpy (cmin, acmin, 3 * sizeof(float));
    memcpy (cmax, acmax, 3 * sizeof(float));

    root = new ADTreeNode3FM;
    for (int i = 0; i < 3; i++)
      root->sep[i] = (cmin[i] + cmax[i]) / 2;
  }

  ADTree3FM :: ~ADTree3FM ()
  {
    root->DeleteChilds();
    delete root;
  }


  void ADTree3FM :: Insert (const float * p, int pi)
  {
    ADTreeNode3FM *node;
    ADTreeNode3FM *next;
    int lr;

    float bmin[3];
    float bmax[3];
    int i, dir;
  
    memcpy (bmin, cmin, 3 * sizeof(float));
    memcpy (bmax, cmax, 3 * sizeof(float));

    next = root;
    while (next)
      {
	node = next;
      
	for (i = 0; i < ADTN_SIZE; i++)
	  if (!node->pi[i])
	    {    
	      memcpy (node->data[i], p, 3 * sizeof(float));
	      node->pi[i] = pi;
	    
	      if (ela.Size() < pi)
		ela.SetSize (pi);
	      ela.Elem(pi) = node;
	    
	      return;
	    }

	dir = 0;
	for (i = 0; i < 3; i++)
	  {
	    if (node->sep[i] > p[i])
	      {
		bmax[i] = node->sep[i];
	      }
	    else
	      {
		bmin[i] = node->sep[i];
		dir += (1 << i);
	      }
	  }
	next = node->childs[dir];

	/*
	  if (node->sep > p[dir])
	  {
	  next = node->left;
	  bmax[dir] = node->sep;
	  lr = 0;
	  }
	  else
	  {
	  next = node->right;
	  bmin[dir] = node->sep;
	  lr = 1;
	  }
	*/
      }


    next = new ADTreeNode3FM;
    memcpy (next->data[0], p, 3 * sizeof(float));
    next->pi[0] = pi;

    for (i = 0; i < 3; i++)
      next->sep[i] = (bmin[i] + bmax[i]) / 2;
  

    if (ela.Size() < pi)
      ela.SetSize (pi);
    ela.Elem(pi) = next;

    node->childs[dir] = next;
    next->father = node;

    while (node)
      {
	node->nchilds++;
	node = node->father;
      }
  }

  void ADTree3FM :: DeleteElement (int pi)
  {
    ADTreeNode3FM * node = ela.Get(pi);

    int i;
    for (i = 0; i < ADTN_SIZE; i++)
      if (node->pi[i] == pi)
	node->pi[i] = 0;

    node = node->father;
    while (node)
      {
	node->nchilds--;
	node = node->father;
      }
  }

  void ADTree3FM :: GetIntersecting (const float * bmin, 
				     const float * bmax,
				     NgArray<int> & pis) const
  {
    static NgArray<ADTreeNode3FM*> stack(1000);
    ADTreeNode3FM * node;
    int dir, i, stacks;

    stack.SetSize (1000);
    pis.SetSize(0);

    stack.Elem(1) = root;
    stacks = 1;

    while (stacks)
      {
	node = stack.Get(stacks);
	stacks--;

	int * hpi = node->pi;
	for (i = 0; i < ADTN_SIZE; i++)
	  if (hpi[i])
	    {
	      float * datai = &node->data[i][0];
	      if (datai[0] >= bmin[0] && datai[0] <= bmax[0] &&
		  datai[1] >= bmin[1] && datai[1] <= bmax[1] &&
		  datai[2] >= bmin[2] && datai[2] <= bmax[2])
	      
		pis.Append (node->pi[i]);
	    }

	/*
	  if (node->pi)
	  {
	  if (node->data[0] >= bmin[0] && node->data[0] <= bmax[0] &&
	  node->data[1] >= bmin[1] && node->data[1] <= bmax[1] &&
	  node->data[2] >= bmin[2] && node->data[2] <= bmax[2])

	  pis.Append (node->pi);
	  }
	*/
      
	int i1min = (bmin[0] <= node->sep[0]) ? 0 : 1;
	int i1max = (bmax[0] < node->sep[0]) ? 0 : 1;
	int i2min = (bmin[1] <= node->sep[1]) ? 0 : 1;
	int i2max = (bmax[1] < node->sep[1]) ? 0 : 1;
	int i3min = (bmin[2] <= node->sep[2]) ? 0 : 1;
	int i3max = (bmax[2] < node->sep[2]) ? 0 : 1;

	int i1, i2, i3;
	for (i1 = i1min; i1 <= i1max; i1++)
	  for (i2 = i2min; i2 <= i2max; i2++)
	    for (i3 = i3min; i3 <= i3max; i3++)
	      {
		i = i1+2*i2+4*i3;
		if (node->childs[i])
		  {
		    stacks++;
		    stack.Elem(stacks) = node->childs[i];
		  }
	      }
      
	/*
	  if (node->left && bmin[dir] <= node->sep)
	  {
	  stacks++;
	  stack.Elem(stacks) = node->left;
	  stackdir.Elem(stacks) = ndir;
	  }
	  if (node->right && bmax[dir] >= node->sep)
	  {
	  stacks++;
	  stack.Elem(stacks) = node->right;
	  stackdir.Elem(stacks) = ndir;
	  }
	*/
      }
  }

  void ADTree3FM :: PrintRec (ostream & ost, const ADTreeNode3FM * node) const
  {
    int i;
    if (node->data)
      {
	ost << node->pi << ": ";
	ost << node->nchilds << " childs, ";
	for (i = 0; i < 3; i++)
	  ost << node->data[i] << " ";
	ost << endl;
      }

    for (i = 0; i < 8; i++)
      if (node->childs[i])
	PrintRec (ost, node->childs[i]);
  }




#endif
  





  /* ******************************* ADTree6 ******************************* */


  ADTreeNode6 :: ADTreeNode6()
  {
    pi = -1;

    left = NULL;
    right = NULL;
    father = NULL;
    nchilds = 0;
  }

  void ADTreeNode6 :: DeleteChilds ()
  {
    if (left)
      {
	left->DeleteChilds();
	delete left;
	left = NULL;
      }
    if (right)
      {
	right->DeleteChilds();
	delete right;
	right = NULL;
      }
  }


  BlockAllocator ADTreeNode6 :: ball (sizeof (ADTreeNode6));
  void * ADTreeNode6 :: operator new(size_t s)
  {
    return ball.Alloc();
  }

  void ADTreeNode6 :: operator delete (void * p)
  {
    ball.Free (p);
  }





  ADTree6 :: ADTree6 (const float * acmin, 
		      const float * acmax)
    : ela(0)
  {
    memcpy (cmin, acmin, 6 * sizeof(float));
    memcpy (cmax, acmax, 6 * sizeof(float));

    root = new ADTreeNode6;
    root->sep = (cmin[0] + cmax[0]) / 2;
  }

  ADTree6 :: ~ADTree6 ()
  {
    root->DeleteChilds();
    delete root;
  }

  void ADTree6 :: Insert (const float * p, int pi)
  {
    ADTreeNode6 *node(NULL);
    ADTreeNode6 *next;
    int dir;
    int lr(0);

    float bmin[6];
    float bmax[6];

  
    memcpy (bmin, cmin, 6 * sizeof(float));
    memcpy (bmax, cmax, 6 * sizeof(float));

    next = root;
    dir = 0;
    while (next)
      {
	node = next;

	if (node->pi == -1)
	  {    
	    memcpy (node->data, p, 6 * sizeof(float));
	    node->pi = pi;

	    if (ela.Size() < pi+1)
	      ela.SetSize (pi+1);
	    ela[pi] = node;

	    return;
	  }

	if (node->sep > p[dir])
	  {
	    next = node->left;
	    bmax[dir] = node->sep;
	    lr = 0;
	  }
	else
	  {
	    next = node->right;
	    bmin[dir] = node->sep;
	    lr = 1;
	  }

	dir++;
	if (dir == 6) dir = 0;
      }


    next = new ADTreeNode6;
    memcpy (next->data, p, 6 * sizeof(float));
    next->pi = pi;
    next->sep = (bmin[dir] + bmax[dir]) / 2;

    if (ela.Size() < pi+1)
      ela.SetSize (pi+1);
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

  void ADTree6 :: DeleteElement (int pi)
  {
    ADTreeNode6 * node = ela[pi];

    node->pi = -1;

    node = node->father;
    while (node)
      {
	node->nchilds--;
	node = node->father;
      }
  }

  void ADTree6 :: PrintMemInfo (ostream & ost) const
  {
    ost << Elements() << " elements a " << sizeof(ADTreeNode6) 
	<< " Bytes = "
	<< Elements() * sizeof(ADTreeNode6) << endl;
    ost << "maxind = " << ela.Size() << " = " << sizeof(ADTreeNode6*) * ela.Size() << " Bytes" << endl;
  }



  class inttn6 {
  public:
    int dir;
    ADTreeNode6 * node;
  };




  void ADTree6 :: GetIntersecting (const float * bmin, 
				   const float * bmax,
				   NgArray<int> & pis) const
  {
    // static NgArray<inttn6> stack(10000);
    // stack.SetSize (10000);
    NgArrayMem<inttn6,10000> stack(10000);
    pis.SetSize(0);

    stack[0].node = root;
    stack[0].dir = 0;
    int stacks = 0;

    while (stacks >= 0)
      {
	ADTreeNode6 * node = stack[stacks].node;
	int dir = stack[stacks].dir; 

	stacks--;
	if (node->pi != -1)
	  {
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
	  }

	int ndir = (dir+1) % 6;

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

  void ADTree6 :: PrintRec (ostream & ost, const ADTreeNode6 * node) const
  {
    
    // if (node->data)     // true anyway
      {
	ost << node->pi << ": ";
	ost << node->nchilds << " childs, ";
	for (int i = 0; i < 6; i++)
	  ost << node->data[i] << " ";
	ost << endl;
      }
    if (node->left)
      PrintRec (ost, node->left);
    if (node->right)
      PrintRec (ost, node->right);
  }


  int ADTree6 :: DepthRec (const ADTreeNode6 * node) const
  {
    int ldepth = 0;
    int rdepth = 0;

    if (node->left)
      ldepth = DepthRec(node->left);
    if (node->right)
      rdepth = DepthRec(node->right);
    return 1 + max2 (ldepth, rdepth);
  }

  int ADTree6 :: ElementsRec (const ADTreeNode6 * node) const
  {
    int els = 1;
    if (node->left)
      els += ElementsRec(node->left);
    if (node->right)
      els += ElementsRec(node->right);
    return els;
  }






  /*
  template <int dim, typename T>
  T_ADTree<dim,T> :: T_ADTree (Point<dim> acmin, Point<dim> acmax)
  // : ela(0)
  {
    cmin = acmin;
    cmax = acmax;
    
    root = new T_ADTreeNode<dim,T>;
    root->sep = (cmin[0] + cmax[0]) / 2;
  }
  */

  /*
  template <int dim, typename T>  
  T_ADTree<dim,T> :: ~T_ADTree ()
  {
    root->DeleteChilds();
    delete root;
  }
  */

  /*
  template <int dim, typename T>
  void T_ADTree<dim,T> :: Insert (Point<dim> p, T pi)
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


    next = new T_ADTreeNode<dim,T>;
    // memcpy (next->data, p, dim * sizeof(float));
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
  */

  /*
  template <int dim, typename T>
  void T_ADTree<dim,T> :: DeleteElement (T pi)
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
  */

  /*
  template <int dim, typename T>
  void T_ADTree<dim,T> :: PrintMemInfo (ostream & ost) const
  {
    ost << Elements() << " elements a " << sizeof(ADTreeNode6) 
	<< " Bytes = "
	<< Elements() * sizeof(T_ADTreeNode<dim,T>) << endl;
    ost << "maxind = " << ela.Size() << " = " << sizeof(T_ADTreeNode<dim,T>*) * ela.Size() << " Bytes" << endl;
  }
  */

  /*
  template <int dim, typename T>
  class inttn {
  public:
    int dir;
    T_ADTreeNode<dim,T> * node;
  };


  template <int dim, typename T>
  void T_ADTree<dim,T> :: GetIntersecting (Point<dim> bmin, Point<dim> bmax,
                                           NgArray<T> & pis) const
  {
    // static NgArray<inttn6> stack(10000);
    // stack.SetSize (10000);
    NgArrayMem<inttn<dim,T>,10000> stack(10000);
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

            //	if (node->data[0] > bmax[0] || 
            // node->data[1] > bmax[1] || 
            // node->data[2] > bmax[2] || 
            // node->data[3] < bmin[3] || 
            // node->data[4] < bmin[4] || 
            // node->data[5] < bmin[5])
            // ;
	    // else
            // {
            // pis.Append (node->pi);
            // }
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
*/

  /*
  template <int dim, typename T>
  void T_ADTree<dim,T> :: PrintRec (ostream & ost, const T_ADTreeNode<dim,T> * node) const
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

  template <int dim, typename T>
  int T_ADTree<dim,T> :: DepthRec (const T_ADTreeNode<dim,T> * node) const
  {
    int ldepth = 0;
    int rdepth = 0;

    if (node->left)
      ldepth = DepthRec(node->left);
    if (node->right)
      rdepth = DepthRec(node->right);
    return 1 + max2 (ldepth, rdepth);
  }

  template <int dim, typename T>
  int T_ADTree<dim,T> :: ElementsRec (const T_ADTreeNode<dim,T> * node) const
  {
    int els = 1;
    if (node->left)
      els += ElementsRec(node->left);
    if (node->right)
      els += ElementsRec(node->right);
    return els;
  }
*/


  



#ifdef ABC

  /* ******************************* ADTree6F ******************************* */


  ADTreeNode6F :: ADTreeNode6F()
  {
    pi = 0;
    father = NULL;
    nchilds = 0;
    int i;
    for (i = 0; i < 64; i++)
      childs[i] = NULL;
  }

  void ADTreeNode6F :: DeleteChilds ()
  {
    int i;

    for (i = 0; i < 64; i++)
      {
	if (childs[i])
	  childs[i]->DeleteChilds();
	delete childs[i];
	childs[i] = NULL;
      }
  }


  BlockAllocator ADTreeNode6F :: ball(sizeof (ADTreeNode6F));

  void * ADTreeNode6F :: operator new(size_t)
  {
    return ball.Alloc();
  }

  void ADTreeNode6F :: operator delete (void * p)
  {
    ball.Free (p);
  }







  ADTree6F :: ADTree6F (const float * acmin, 
			const float * acmax)
    : ela(0)
  {
    memcpy (cmin, acmin, 6 * sizeof(float));
    memcpy (cmax, acmax, 6 * sizeof(float));

    root = new ADTreeNode6F;
    for (int i = 0; i < 6; i++)
      root->sep[i] = (cmin[i] + cmax[i]) / 2;
  }

  ADTree6F :: ~ADTree6F ()
  {
    root->DeleteChilds();
    delete root;
  }


  void ADTree6F :: Insert (const float * p, int pi)
  {
    ADTreeNode6F *node;
    ADTreeNode6F *next;
    int lr;

    float bmin[6];
    float bmax[6];
    int i, dir;
  
    memcpy (bmin, cmin, 6 * sizeof(float));
    memcpy (bmax, cmax, 6 * sizeof(float));

    next = root;
    while (next)
      {
	node = next;
      
	if (!node->pi)
	  {    
	    memcpy (node->data, p, 6 * sizeof(float));
	    node->pi = pi;

	    if (ela.Size() < pi)
	      ela.SetSize (pi);
	    ela.Elem(pi) = node;

	    return;
	  }

	dir = 0;
	for (i = 0; i < 6; i++)
	  {
	    if (node->sep[i] > p[i])
	      {
		bmax[i] = node->sep[i];
	      }
	    else
	      {
		bmin[i] = node->sep[i];
		dir += (1 << i);
	      }
	  }
	next = node->childs[dir];

	/*
	  if (node->sep > p[dir])
	  {
	  next = node->left;
	  bmax[dir] = node->sep;
	  lr = 0;
	  }
	  else
	  {
	  next = node->right;
	  bmin[dir] = node->sep;
	  lr = 1;
	  }
	*/
      }


    next = new ADTreeNode6F;
    memcpy (next->data, p, 6 * sizeof(float));
    next->pi = pi;

    for (i = 0; i < 6; i++)
      next->sep[i] = (bmin[i] + bmax[i]) / 2;
  

    if (ela.Size() < pi)
      ela.SetSize (pi);
    ela.Elem(pi) = next;

    node->childs[dir] = next;
    next->father = node;

    while (node)
      {
	node->nchilds++;
	node = node->father;
      }
  }

  void ADTree6F :: DeleteElement (int pi)
  {
    ADTreeNode6F * node = ela.Get(pi);

    node->pi = 0;

    node = node->father;
    while (node)
      {
	node->nchilds--;
	node = node->father;
      }
  }

  void ADTree6F :: GetIntersecting (const float * bmin, 
				    const float * bmax,
				    NgArray<int> & pis) const
  {
    static NgArray<ADTreeNode6F*> stack(1000);
    ADTreeNode6F * node;
    int dir, i, stacks;

    stack.SetSize (1000);
    pis.SetSize(0);

    stack.Elem(1) = root;
    stacks = 1;

    while (stacks)
      {
	node = stack.Get(stacks);
	stacks--;

	if (node->pi)
	  {
	    if (
		node->data[0] >= bmin[0] && node->data[0] <= bmax[0] &&
		node->data[1] >= bmin[1] && node->data[1] <= bmax[1] &&
		node->data[2] >= bmin[2] && node->data[2] <= bmax[2] &&
		node->data[3] >= bmin[3] && node->data[3] <= bmax[3] &&
		node->data[4] >= bmin[4] && node->data[4] <= bmax[4] &&
		node->data[5] >= bmin[5] && node->data[5] <= bmax[5]
		)

	      pis.Append (node->pi);
	  }

      
	int i1min = (bmin[0] <= node->sep[0]) ? 0 : 1;
	int i1max = (bmax[0] < node->sep[0]) ? 0 : 1;
	int i2min = (bmin[1] <= node->sep[1]) ? 0 : 1;
	int i2max = (bmax[1] < node->sep[1]) ? 0 : 1;
	int i3min = (bmin[2] <= node->sep[2]) ? 0 : 1;
	int i3max = (bmax[2] < node->sep[2]) ? 0 : 1;

	int i4min = (bmin[3] <= node->sep[3]) ? 0 : 1;
	int i4max = (bmax[3] <  node->sep[3]) ? 0 : 1;
	int i5min = (bmin[4] <= node->sep[4]) ? 0 : 1;
	int i5max = (bmax[4] <  node->sep[4]) ? 0 : 1;
	int i6min = (bmin[5] <= node->sep[5]) ? 0 : 1;
	int i6max = (bmax[5] <  node->sep[5]) ? 0 : 1;

	int i1, i2, i3, i4, i5, i6;
	for (i1 = i1min; i1 <= i1max; i1++)
	  for (i2 = i2min; i2 <= i2max; i2++)
	    for (i3 = i3min; i3 <= i3max; i3++)
	      for (i4 = i4min; i4 <= i4max; i4++)
		for (i5 = i5min; i5 <= i5max; i5++)
		  for (i6 = i6min; i6 <= i6max; i6++)
		    {
		      i = i1 + 2*i2 + 4*i3 + 8*i4 + 16*i5 +32*i6;
		      if (node->childs[i])
			{
			  stacks++;
			  stack.Elem(stacks) = node->childs[i];
			}
		    }
      
	/*
	  if (node->left && bmin[dir] <= node->sep)
	  {
	  stacks++;
	  stack.Elem(stacks) = node->left;
	  stackdir.Elem(stacks) = ndir;
	  }
	  if (node->right && bmax[dir] >= node->sep)
	  {
	  stacks++;
	  stack.Elem(stacks) = node->right;
	  stackdir.Elem(stacks) = ndir;
	  }
	*/
      }
  }

  void ADTree6F :: PrintRec (ostream & ost, const ADTreeNode6F * node) const
  {
    int i;
    if (node->data)
      {
	ost << node->pi << ": ";
	ost << node->nchilds << " childs, ";
	for (i = 0; i < 6; i++)
	  ost << node->data[i] << " ";
	ost << endl;
      }

    for (i = 0; i < 64; i++)
      if (node->childs[i])
	PrintRec (ost, node->childs[i]);
  }



#endif



  /* ************************************* Point3dTree ********************** */



  Point3dTree :: Point3dTree (const Point<3> & pmin, const Point<3> & pmax)
  {
    float pmi[3], pma[3];
    for (int i = 0; i < 3; i++)
      {
	pmi[i] = pmin(i);
	pma[i] = pmax(i);
      }
    tree = new ADTree3 (pmi, pma);
  }

  Point3dTree :: ~Point3dTree ()
  {
    delete tree;
  }



  void Point3dTree :: Insert (const Point<3> & p, int pi)
  {
    float pd[3];
    pd[0] = p(0);
    pd[1] = p(1);
    pd[2] = p(2);
    tree->Insert (pd, pi);
  }

  void Point3dTree :: GetIntersecting (const Point<3> & pmin, const Point<3> & pmax, 
				       NgArray<int> & pis) const
  {
    float pmi[3], pma[3];
    for (int i = 0; i < 3; i++)
      {
	pmi[i] = pmin(i);
	pma[i] = pmax(i);
      }
    tree->GetIntersecting (pmi, pma, pis);
  }







  /*
  template <int dim, typename T>
  BoxTree<dim,T> :: BoxTree (const Box<dim> & abox)
  {
    boxpmin = abox.PMin();
    boxpmax = abox.PMax();
    Point<2*dim> tpmin, tpmax;
    for (int i = 0; i < dim; i++)
      {
	tpmin(i) = tpmin(i+dim) = boxpmin(i);
	tpmax(i) = tpmax(i+dim) = boxpmax(i);
      }
    tree = new T_ADTree<2*dim,T> (tpmin, tpmax);
  }
  */

  /*
  template <int dim, typename T>
  BoxTree<dim,T> :: BoxTree (const Point<dim> & apmin, const Point<dim> & apmax)
  {
    boxpmin = apmin;
    boxpmax = apmax;
    Point<2*dim> tpmin, tpmax;
    for (int i = 0; i < dim; i++)
      {
	tpmin(i) = tpmin(i+dim) = boxpmin(i);
	tpmax(i) = tpmax(i+dim) = boxpmax(i);
      }
    tree = new T_ADTree<2*dim,T> (tpmin, tpmax);
  }

  template <int dim, typename T>
  BoxTree<dim,T> :: ~BoxTree ()
  {
    delete tree;
  }

  template <int dim, typename T>
  void BoxTree<dim,T> :: Insert (const Point<dim> & bmin, const Point<dim> & bmax, T pi)
  {
    Point<2*dim> tp;

    for (size_t i = 0; i < dim; i++)
      {
	tp(i) = bmin(i);
	tp(i+dim) = bmax(i);
      }

    tree->Insert (tp, pi);
  }
  */

  /*
  template <int dim, typename T>
  void BoxTree<dim,T> ::GetIntersecting (const Point<dim> & pmin, const Point<dim> & pmax, 
                                         NgArray<T> & pis) const
  {
    Point<2*dim> tpmin, tpmax;
    double tol = Tolerance();
    for (size_t i = 0; i < dim; i++)
      {
	tpmin(i) = boxpmin(i);
	tpmax(i) = pmax(i)+tol;
        
	tpmin(i+dim) = pmin(i)-tol;
	tpmax(i+dim) = boxpmax(i);
      }

    tree->GetIntersecting (tpmin, tpmax, pis);
  }
  */
  
  /*
  template<> BlockAllocator T_ADTreeNode<4,INDEX> :: ball(sizeof (T_ADTreeNode<4,INDEX>));
  template class T_ADTree<4,INDEX>;
  template class BoxTree<2,INDEX>;

  template<> BlockAllocator T_ADTreeNode<4,INDEX_2> :: ball(sizeof (T_ADTreeNode<4,INDEX_2>));
  template class T_ADTree<4,INDEX_2>;
  template class BoxTree<2,INDEX_2>;
  
  
  template<> BlockAllocator T_ADTreeNode<6,INDEX> :: ball(sizeof (T_ADTreeNode<6,INDEX>));
  template class T_ADTree<6,INDEX>;
  template class BoxTree<3,INDEX>;
  */
  
}

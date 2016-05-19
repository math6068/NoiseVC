/************************************************
** This is a local search solver for Minimum Vertex Cover.                                                       
************************************************/


/************************************************
** Date:	2016.4.20  
** NoiseVC     
** Author: Zongjie Ma, Yi Fan, Chengqian Li                                                                     
************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <sys/times.h>
#include <cmath>

using namespace std;

#define fastvc_mode 0
#define heap_cover_mode 1
#define one_tabu_remove_mode 0
#define one_tabu_add_mode 0
#define add_time_stamp_mode 0
#define remove_time_stamp_mode 0

//#define debug_mode

#define pop(stack) stack[--stack ## _fill_pointer]
#define push(item, stack) stack[stack ## _fill_pointer++] = item

/*max vertex count and max edge count*/
#define	MAXV	15000000
#define MAXE	80000000

#define individual_analysis_on_init_sls_mode

#ifdef individual_analysis_on_init_sls_mode
double	init_time;
double	sls_time;
double	sls_step_speed_per_ms;
#endif

tms start, finish;
int start_time;

struct Edge{
	int v1;
	int v2;
};

/*parameters of algorithm*/
long long	max_steps;			//step limit
int			cutoff_time;			//time limit
long long	step;
int			optimal_size;			//terminate the algorithm before step limit if it finds a vertex cover of optimal_size

/*parameters of the instance*/
int		v_num;//|V|: 1...v
int		e_num;//|E|: 0...e-1

/*structures about edge*/
Edge	edge[MAXE];  

/*structures about vertex*/
int		dscore[MAXV];			//dscore of v
long long	time_stamp[MAXV];


//from vertex to it's edges and neighbors
int*	v_edges[MAXV];	//edges related to v, v_edges[i][k] means vertex v_i's k_th edge
int*	v_adj[MAXV];		//v_adj[v_i][k] = v_j(actually, that is v_i's k_th neighbor)
int		v_degree[MAXV];	//amount of edges (neighbors) related to v


/* structures about solution */
//current candidate solution
int		c_size;						//cardinality of C
bool	v_in_c[MAXV];				//a flag indicates whether a vertex is in C
#if 0
int		remove_cand[MAXV];			//remove candidates, an array consists of only vertices in C, not including tabu_remove


int		index_in_remove_cand[MAXV];
int		remove_cand_size;
#endif

#if heap_cover_mode == 1
int		heap_vc_size;
int		heap_vc[MAXV];
int		index_in_heap_vc[MAXV];
#endif

//best solution found
int		best_c_size;
bool	best_v_in_c[MAXV];			//a flag indicates whether a vertex is in best solution
double  best_comp_time;
long    best_step;


//uncovered edge stack
int		uncov_stack[MAXE];		//store the uncov edge number
int		uncov_stack_fill_pointer;
int		index_in_uncov_stack[MAXE];//which position is an edge in the uncov_stack


//CC and taboo
//int 	conf_change[MAXV];
//int		tabu_remove=0;
#if one_tabu_remove_mode == 1
int tabu_remove_v;
#endif

#if one_tabu_add_mode == 1
int	tabu_add_v;
#endif

#if add_time_stamp_mode == 1
long long	add_stamp[MAXV];
#endif

#if remove_time_stamp_mode == 1
long long	remove_stamp[MAXV];
#endif


/* functions declaration */
int build_instance(char *filename);
void init_sol();
void cover_LS();
void add(int v);
void remove(int v);
#if heap_cover_mode == 1
void Error(char* s);
bool heap_worse(const int &a, const int &b);
void heap_up(int p);
void heap_down(int p);
void heap_add(int v);
void heap_del(int v);
int checkHeapCoverMap();
#endif
void update_edge_weight();
void cover_rest_edges();
int check_solution();



void update_best_sol()
{
	int i;
#if fastvc_mode == 1
	for (i=1;i<=v_num;i++)
		best_v_in_c[i] = v_in_c[i];
	best_c_size = c_size;
#endif

#if heap_cover_mode == 1
	int v;
	for(v = 1; v <= v_num; v++)
	{
		if(index_in_heap_vc[v])
			best_v_in_c[v] = 1;
		else
			best_v_in_c[v] = 0;
	}
	best_c_size = heap_vc_size;
#endif	
	
	times(&finish);
	best_comp_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
	best_comp_time = round(best_comp_time * 100)/100.0;
	best_step = step;
}



int v_degree_tmp[MAXV];

int build_instance(char *filename)
{
	char line[1024];
	char tempstr1[10];
	char tempstr2[10];
	int  v,e;
	
	char	tmp;
	int		v1,v2;
	
	ifstream infile(filename);
    if(infile==NULL) return 0;

	/*** build problem data structures of the instance ***/
	infile.getline(line,1024);
	while (line[0] != 'p') infile.getline(line,1024);
	sscanf(line, "%s %s %d %d", tempstr1, tempstr2, &v_num, &e_num);
	
	if (v_num > MAXV) {
		cout<<"the number of vertices ("<<v_num<<") exceeds MAXV ("<<MAXV<<")."<<endl;
		exit(0);
	}
	if (e_num > MAXE) {
		cout<<"the number of vertices ("<<e_num<<") exceeds MAXV ("<<MAXE<<")."<<endl;
		exit(0);
	}

	/* read edges and compute v_degree */
	for (v=1; v<=v_num; v++) v_degree[v] = 0;
	
	for (e=0; e<e_num; e++)
	{
		infile>>tmp>>v1>>v2;
		v_degree[v1]++;
		v_degree[v2]++;
		
		edge[e].v1 = v1;
		edge[e].v2 = v2;
	}
	infile.close();
	
	/* build v_adj and v_edges arrays */
	for (v=1; v<=v_num; v++)
	{
		v_adj[v] = new int[v_degree[v]];
		v_edges[v] = new int[v_degree[v]];
	}
	
	//for(v=1; v<=v_num; v++) v_degree_tmp[v]=0;
	
	for (e=0; e<e_num; e++)
	{
		v1=edge[e].v1;
		v2=edge[e].v2;

		v_edges[v1][v_degree_tmp[v1]] = e;
		v_edges[v2][v_degree_tmp[v2]] = e;

		v_adj[v1][v_degree_tmp[v1]] = v2;
		v_adj[v2][v_degree_tmp[v2]] = v1;

		v_degree_tmp[v1]++;
		v_degree_tmp[v2]++;
	}
	
	return 1;
}


void free_memory()
{
	for (int v=1; v<=v_num; v++)
	{
		delete[] v_adj[v];
		delete[] v_edges[v];
	}
}

void reset_remove_cand()
{
/*
	int v,j;
	j=0;
	for (v=1; v<=v_num; v++)
	{
		if(v_in_c[v]==1)
		{
			remove_cand[j] = v;
			index_in_remove_cand[v]=j;
			j++;
		}
		else index_in_remove_cand[v]=0;
	}
	
	remove_cand_size = j;
*/
}



void update_target_size()
{	
	int v,i;
	int best_dscore;
	int best_remove_v;//vertex with the highest improvement in C
#if heap_cover_mode == 1
	best_remove_v = heap_vc[1];
	heap_del(best_remove_v);
#endif
#if fastvc_mode == 1
	c_size--;
	best_remove_v = remove_cand[0];
	best_dscore = dscore[best_remove_v];
	
	if(dscore[best_remove_v]!=0)
	{
		for (i=1; i<remove_cand_size; ++i)
		{
			v = remove_cand[i];
			
			if(dscore[v]==0) break;

			if (dscore[v] > dscore[best_remove_v])
				best_remove_v = v;
		}
	}
#endif
/*#if heap_cover_mode == 1
	//int v,i;
	//int best_dscore; kaolvjiashang
	//int best_remove_v;//vertex with the highest improvement in C
	//best_remove_v = heap_vc[1];
	//heap_del(best_remove_v);
int cand_counts=50;
	int i,v;

	int best_remove_v = heap_vc[rand()%heap_vc_size+1];
	
	for (i=1; i<cand_counts; ++i)
	{
		v = heap_vc[rand()%heap_vc_size+1];
	
		if( dscore[v] < dscore[best_remove_v])
			continue;
		else if( dscore[v]> dscore[best_remove_v] )
			best_remove_v = v;
		else if (time_stamp[v]<time_stamp[best_remove_v])
			best_remove_v = v;
	} 
     heap_del(best_remove_v);
#endif	*/
	remove(best_remove_v);
#if fastvc_mode == 1	
	//remove best_remove_v from remove_cand, and move the last vertex in remove_cand to the position
#if 0
	int last_remove_cand_v = remove_cand[--remove_cand_size];
	int index = index_in_remove_cand[best_remove_v];
	remove_cand[index] = last_remove_cand_v;
	index_in_remove_cand[last_remove_cand_v] = index;
#endif
#endif
	//reset_remove_cand();
}




//update the best vertex in C 
int cand_count = 50;


int choose_remove_v()
{
#if fastvc_mode == 1
	int i,v;

		int best_v = remove_cand[rand()%remove_cand_size];
	
		for (i=1; i<cand_count; ++i)
		{
			v = remove_cand[rand()%remove_cand_size];
		
			if( dscore[v] < dscore[best_v])
				continue;
			else if( dscore[v]> dscore[best_v] )
				best_v = v;
			else if (time_stamp[v]<time_stamp[best_v])
				best_v = v;
		}
		
		return best_v;
#endif

#if heap_cover_mode == 1     //along with the next 2 line by me 

/*if (dscore[heap_vc[1]] == 0)   //del 0 mode
    return heap_vc[1];
else if (rand() % 1000 < 0 )
	return heap_vc[1];
else
    return heap_vc[rand()%heap_vc_size+1];*/
if (rand() % 1000 < k )      // k is set to 400 in this work
	return heap_vc[1];
else
    return heap_vc[rand()%heap_vc_size+1];
    
#endif

/*#if heap_cover_mode == 1
	#if one_tabu_remove_mode == 1
	if(heap_vc[1] != tabu_remove_v) return heap_vc[1];
#ifdef debug_mode
	cout << "tabu functions" << endl;
#endif
	if(heap_worse(heap_vc[2], heap_vc[3])) return heap_vc[3];
	return heap_vc[2];
	#else
	return heap_vc[1];
	#endif
#endif*/
}



inline
void uncover(int e) 
{
	index_in_uncov_stack[e] = uncov_stack_fill_pointer;
	push(e,uncov_stack);
}


inline
void cover(int e)
{
	int index,last_uncov_edge;

	//since the edge is satisfied, its position can be reused to store the last_uncov_edge
	last_uncov_edge = pop(uncov_stack);
	index = index_in_uncov_stack[e];
	uncov_stack[index] = last_uncov_edge;
	index_in_uncov_stack[last_uncov_edge] = index;
}



//int uncov_degree[MAXV];

void init_sol()
{
	int i,v,e;
	int v1, v2;

	/*** build solution data structures of the instance ***/
	for (v=1; v<=v_num; v++)
	{
	
		//cout<<"d "<<v_degree[v]<<endl;
#if fastvc_mode == 1
		v_in_c[v] = 0;
#endif
		dscore[v] = 0;
		//conf_change[v] = 1;
		time_stamp[v]= 0; 

#if add_time_stamp_mode == 1
		add_stamp[v] = 0;
#endif

#if remove_time_stamp_mode == 1
		remove_stamp[v] = 0;
#endif
		
		//uncov_degree[v] = v_degree[v];
	}
	
	//construct a vertex cover
#if fastvc_mode == 1
	c_size = 0;
#endif
// init an empty heap
#if heap_cover_mode == 1
	heap_vc_size = 0;
	for(v = 1; v <= v_num; v++)
	{
		index_in_heap_vc[v] = 0;
	}
#endif
	for (e=0; e<e_num; e++)
	{	
		v1=edge[e].v1;
		v2=edge[e].v2;
#if fastvc_mode == 1		
		if (v_in_c[v1]==0 && v_in_c[v2]==0)//if uncovered, choose the endpoint with higher degree
#endif
#if heap_cover_mode == 1
		if(index_in_heap_vc[v1]==0 && index_in_heap_vc[v2]==0)
#endif
		{
#if fastvc_mode == 1
			if(v_degree[v1] > v_degree[v2]) 
			{
				v_in_c[v1]=1;
			}
			else{
				v_in_c[v2]=1;
			}
			c_size++;
#endif
#if heap_cover_mode == 1
			if(v_degree[v1] > v_degree[v2]) 
			{
				v = v1;
			}
			else{
				v = v2;
			}
			heap_vc[++heap_vc_size] = v;
			index_in_heap_vc[v] = heap_vc_size;
#endif		
		}
	}

	//init uncovered edge stack
	uncov_stack_fill_pointer = 0;
	
	//calculate dscores
	for (e=0; e<e_num; e++)
	{
		v1=edge[e].v1;
		v2=edge[e].v2;
#if fastvc_mode == 1		
		if (v_in_c[v1]==1 && v_in_c[v2]==0) dscore[v1]--;
		else if (v_in_c[v2]==1 && v_in_c[v1]==0) dscore[v2]--;
#endif
#if heap_cover_mode == 1
		if (index_in_heap_vc[v2] == 0) dscore[v1]--;
		else if (index_in_heap_vc[v1] == 0) dscore[v2]--;
#endif
	}
	

	//remove redundent vertices
#if fastvc_mode == 1
		for (v=1; v<=v_num; v++)
		{
			if (v_in_c[v]==1 && dscore[v]==0) 
			{
				remove(v);
				c_size--;
			}
		}
#endif
#if heap_cover_mode == 1
	for(i = 1; i <= heap_vc_size; i++)
	{
		v = heap_vc[i];
		if(dscore[v] == 0)
		{
			int last_vertex_in_heap_vc = heap_vc[heap_vc_size--];
			heap_vc[i] = last_vertex_in_heap_vc;
			index_in_heap_vc[last_vertex_in_heap_vc] = i;
			index_in_heap_vc[v] = 0;//should be 0
			i--;//should exist
			for(int j = 0; j < v_degree[v]; j++)
			{
				dscore[v_adj[v][j]]--;
			}
		}
	}
	// establish the heap
	for(i = heap_vc_size; i > 0; --i)
		heap_down(i);
#endif
	
	update_best_sol();//initialize the best found solution

#if fastvc_mode == 1
	//cout<<c_size<<' '<<best_comp_time<<' ';
#if 0
	reset_remove_cand();
#endif
#endif
}


void add(int v)
{
#if fastvc_mode == 1
	v_in_c[v] = 1;
#endif
	dscore[v] = -dscore[v];
	
	int i,e,n;

	int edge_count = v_degree[v];
	
	for (i=0; i<edge_count; ++i)
	{
		e = v_edges[v][i];// v's i'th edge
		n = v_adj[v][i];//v's i'th neighbor
#if fastvc_mode == 1
		if (v_in_c[n]==0)//this adj isn't in cover set
#endif
#if heap_cover_mode == 1
		if(index_in_heap_vc[n] == 0)
#endif
		{
			dscore[n]--;
			//conf_change[n] = 1;

			cover(e);
		}
		else
		{
			dscore[n]++; 
#if heap_cover_mode == 1
			heap_up(index_in_heap_vc[n]);
#endif
		}
	}
	
}

void remove(int v)
{
#if fastvc_mode == 1
	v_in_c[v] = 0;
#endif
	dscore[v] = -dscore[v];
	//conf_change[v] = 0;

	int i,e,n;

	int edge_count = v_degree[v];
	for (i=0; i<edge_count; ++i)
	{
		e = v_edges[v][i];
		n = v_adj[v][i];
#if fastvc_mode == 1
		if (v_in_c[n]==0)//this adj isn't in cover set
#endif
#if heap_cover_mode == 1
		if(index_in_heap_vc[n] == 0)
#endif
		{
			dscore[n]++;  
			//conf_change[n] = 1;

			uncover(e);
		}
		else
		{
			dscore[n]--; 
#if heap_cover_mode == 1
			heap_down(index_in_heap_vc[n]);
#endif
		}
	}
}


#if heap_cover_mode
inline bool heap_worse(const int &a, const int &b)
// if a is worse than b, return 1, otherwise return 0
{
/*
#if add_time_stamp_mode == 1// the same as the common age
	if(dscore[a] < dscore[b] || (dscore[a] == dscore[b] && add_stamp[a] > add_stamp[b]))
*/
#if remove_time_stamp_mode == 1
	if(dscore[a] < dscore[b] || (dscore[a] == dscore[b] && remove_stamp[a] > remove_stamp[b]))
#else
	if(dscore[a] < dscore[b] || (dscore[a] == dscore[b] && time_stamp[a] > time_stamp[b])) 
#endif
		return 1;
	return 0;
}

inline void heap_up(int p)
//used when p may possibly be better than its parent
//push up the item pointed by p to the proper location
//if pushing up operations are not needed, it can return immediately without doing anything
{
#ifdef debug_mode
	if(p <= 0 || p > heap_vc_size)
		Error("gjj308yuuh");
#endif
	int q = p >> 1;// q is the parent node, p is one of the child node
	int a = heap_vc[p];// to be placed in an upper location
	int z = heap_vc[q];// 'z' is above 'a'
	while(q)// still a location in the heap
	{
		if(heap_worse(z, a))// 'z' is worse than 'a'
		{
			// push down 'z'
			heap_vc[p] = z;
			index_in_heap_vc[z] = p;
		}
		else
			break;
		// move 'p' and 'q' up to point to upper locations
		p = q;
		q = p >> 1; // shift one bit position.
		// obtain a new item for next comparisons
		z = heap_vc[q];
	}
	// find the proper location for 'a'
	heap_vc[p] = a;
	index_in_heap_vc[a] = p;
}

inline void heap_down(int p)
//used when p may possibly be worse than either child
//push down the item pointed by p the the proper location
//if pushing down operations are not needed, it can return immediately without doing anything
{
#ifdef debug_mode
	if(p <= 0 || p > heap_vc_size)
		Error("fejgj483yuhj");
#endif

	int q = p << 1;// p is the parent node, q is one of the child nodes 
	int a = heap_vc[p];// parent item
	int x = heap_vc[q];// left child
	int y = heap_vc[q+1];// right child
	
	while(q <= heap_vc_size)// q is a location in the heap
	{
		if(q < heap_vc_size && heap_worse(x, y))// left child is worse, so only the right child can possibly be pushed up
		{
			if(heap_worse(a, y))// the right child is better than its parent
			{
				// push up 'y'
				heap_vc[p] = y;
				index_in_heap_vc[y] = p;
			}
			else
				break;
			// move 'p' and 'q' to point to lower locations
			p = q + 1;
			q = p << 1;		
			// obtain items for next comparisons
			x = heap_vc[q];
			y = heap_vc[q+1];
		}		
		else
		{
			if(heap_worse(a, x))// the right child is worse
			{
				// push up 'x'
				heap_vc[p] = x;
				index_in_heap_vc[x] = p;
			}
			else
				break;
			// move 'p' and 'q' to point to lower locations
			p = q;
			q = p << 1;	
			// obtain items for next comparisons	
			x = heap_vc[q];
			y = heap_vc[q+1];
		}
	}
	heap_vc[p] = a;
	index_in_heap_vc[a] = p;
}

inline void heap_del(int name)
{
//	printf("(((((((((((((((( heap_del %d %d\n", name, heap_vc[1]);
#ifdef debug_mode
	if(!index_in_heap_vc[name])
		Error("jgt49ghj30uyh");
#endif

	if(heap_vc_size == index_in_heap_vc[name])// when deleting the end element
	{
		heap_vc_size--;
		index_in_heap_vc[name] = 0;
		return;
	}

//printf("((((((((((((((((3333333 heap_del %d %d\n",name,heap_vc[1]);
		
	int p = index_in_heap_vc[name];
	
	// erasing the element
	index_in_heap_vc[name] = 0;
//	printf("(((((((((((((((( heap_del %d %d\n",name,index_in_heap_vc[name]);
	heap_vc[p] = heap_vc[heap_vc_size--];
	index_in_heap_vc[heap_vc[p]] = p;		

//	printf("(((((((((((((((( heap_del %d %d\n",name,index_in_heap_vc[name]);

	// keep the structures of the heap
	// each item in the heap can possibly be deleted, so there are the following two cases:

	// if the previous end item is worse than one of its children
	heap_down(p);
	// if the prevous end item is better than its parent
	heap_up(p);
	// only one of the two above statements will be executed
	
//	printf("(((((((((((((((( heap_del %d %d\n",name,index_in_heap_vc[name]);
}

inline void heap_add(int name)
{
	heap_vc[++heap_vc_size] = name;
	index_in_heap_vc[name] = heap_vc_size;
	heap_up(heap_vc_size);
}


void Error(char s[])
{
	cout << "Error[" << s << "]" << endl;
}
#endif


/*On solution*/

void print_solution()
{
	int mis_vertex_count=0;
	
	for (int i=1; i<=v_num; i++)
	{
		if (best_v_in_c[i]!=1)
			mis_vertex_count++;
	}
	
	if(mis_vertex_count+best_c_size!=v_num)
		cout<<"The size of independent set + the size of vertex cover is not equal to |V(G)|!"<<endl;
	
	cout<<"c Best found independent set size = "<<mis_vertex_count<<endl;
	cout<<"c The following output is the found independent set."<<endl;


	for (int i=1; i<=v_num; i++)
	{
		if (best_v_in_c[i]!=1)//output max independent set
			cout<<i<<'\t';
	}
	cout<<endl;

}

//check whether the solution found is a proper solution
int check_solution()
{
	int v, e;
	int actual_c_size = 0;
//check cover-size
	for(v = 1; v <= v_num; ++v)
		if(best_v_in_c[v]) actual_c_size++;
	if(actual_c_size != best_c_size)
	{
		cout << "best_c_size incorrect" << endl;
		return 0;
	}
//check edges
	for(e=0; e<e_num; ++e)
	{
		if(best_v_in_c[edge[e].v1]!=1 && best_v_in_c[edge[e].v2]!=1)
		{
			cout<<"uncovered edge "<<e<<endl;
			return 0;
		}
	}
	
	return 1;
}


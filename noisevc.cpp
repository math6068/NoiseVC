#include "noisevc.h"

int try_step=100;


int edge_cand;

void cover_LS()
{
	int		remove_v, add_v;
	int		remove_dscore, add_dscore;
	int		e,v1,v2;
	int		i;

	step = 1;

#if one_tabu_remove_mode == 1
	tabu_remove_v = 0;
#endif

#if one_tabu_add_mode == 1
	tabu_add_v = 0;
#endif

	while(1)
	{
		if (uncov_stack_fill_pointer == 0)//update best solution if needed
		{
			update_best_sol();
			
			//if (c_size==optimal_size) return;
				
			update_target_size();
			
			continue;
		}
		
		if(step%try_step==0) //check cutoff
		{
			times(&finish);
			double elap_time = (finish.tms_utime + finish.tms_stime - start_time)/sysconf(_SC_CLK_TCK);
			if(elap_time >= cutoff_time)return;
		}
		
		remove_v = choose_remove_v();
		//remove_dscore = dscore[remove_v];
		
		remove(remove_v);
#if heap_cover_mode == 1
		heap_del(remove_v);
#endif
		
		e = uncov_stack[rand()%uncov_stack_fill_pointer];
		v1 = edge[e].v1;
		v2 = edge[e].v2;
#if one_tabu_add_mode == 1                                   //********** wodexiugai add" ==1
			if(v1 == tabu_add_v) 
			{
#ifdef debug_mode
cout << "tabu functions" << endl;
#endif
				add_v = v2;
			}
			else if(v2 == tabu_add_v) 
			{
#ifdef debug_mode
cout << "tabu functions" << endl;
#endif
				add_v = v1;
			}
			else 
#endif	
#if add_time_stamp_mode == 1	
			if(dscore[v1]>dscore[v2] || (dscore[v1]==dscore[v2] && add_stamp[v1]<add_stamp[v2]) )
#else
			if(dscore[v1]>dscore[v2] || (dscore[v1]==dscore[v2] && time_stamp[v1]<time_stamp[v2]) )
#endif
				add_v=v1;
			else add_v=v2;


		add(add_v);
		
#if 0		
		int index = index_in_remove_cand[remove_v];
		index_in_remove_cand[remove_v] = 0;
		
		remove_cand[index] = add_v;
		index_in_remove_cand[add_v] = index;
#endif
		
		time_stamp[add_v]=time_stamp[remove_v]=step;
#if add_time_stamp_mode == 1
		add_stamp[add_v] = step;
#endif
#if remove_time_stamp_mode == 1
		remove_stamp[remove_v] = step;
#endif

// danger: heap_add should never be executed before updating the time stamps
#if heap_cover_mode == 1
		heap_add(add_v);
#endif
#if one_tabu_remove_mode == 1
	tabu_remove_v = add_v;
#endif
#if one_tabu_add_mode == 1
	tabu_add_v = remove_v;
#endif
		//tabu_remove = add_v;
		
		//update_edge_weight();
		
		step++;
	}

}

int main(int argc, char* argv[])
{
	int seed,i;

	//cout<<"c This is NuMVC, a local search solver for the Minimum Vertex Cover (and also Maximum Independent Set) problem."<<endl;
	
	if(build_instance(argv[1])!=1){
		cout<<"can't open instance file"<<endl;
		return -1;
	}
		optimal_size=0;
		i=2;
		//sscanf(argv[i++],"%d",&cand_count);//if you want to stop the algorithm only cutoff time is reached, set optimal_size to 0.
		//sscanf(argv[i++],"%d",&edge_cand);
		sscanf(argv[i++],"%d",&seed);
		sscanf(argv[i++],"%d",&cutoff_time);

	
		srand(seed);

		//cout<<seed<<' ';
		//cout<<argv[1]<<' ';
		
		times(&start);
		start_time = start.tms_utime + start.tms_stime;

    	init_sol();

#ifdef individual_analysis_on_init_sls_mode
		times(&finish);
		init_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime) / sysconf(_SC_CLK_TCK);
		init_time = round(init_time * 100)/100.0;
#endif

		//if(c_size + uncov_stack_fill_pointer > optimal_size ) 
		//{
			//cout<<"c Start local search..."<<endl;
			cover_LS();
		//}
#ifdef individual_analysis_on_init_sls_mode
		times(&finish);
		sls_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime) / sysconf(_SC_CLK_TCK) - init_time;
		sls_step_speed_per_ms = double(step) * 0.001 / sls_time;
#endif			
		//check solution
		if(check_solution()==1)
		{
			cout << "o " << best_c_size << endl;
			//cout << best_c_size << ' ';

			//print_mvc_solution();
			cout << "c searchSteps " << best_step << endl;
			//printf("%ld ", best_step);
			cout << "c solveTime " << best_comp_time << endl;
			//cout << "c stepVelocity(/0.001ms) " << (long double)(best_step) / (best_comp_time * 1000000) << endl;
			/*cout<<"c Best found vertex cover size = "<<best_c_size<<endl;
			print_solution();
			cout<<"c searchSteps = "<<best_step<<endl;
			cout<<"c solveTime = "<<best_comp_time<<endl;*/
			
			//cout<<best_c_size<<' '<<best_comp_time<<' '<<best_step<<endl;
#ifdef 	individual_analysis_on_init_sls_mode
		//cout<<"c initTime " << init_time << endl;
		//cout<<"c slsTime " << sls_time << endl;
		cout<<"c stepSpeed(/ms) "<< sls_step_speed_per_ms << endl;
#endif
		}
		else
		{
			cout<<"the solution is wrong."<<endl;
			//print_solution();
		}
	
		free_memory();

	return 0;
}

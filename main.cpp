#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "SimStruct.h"
#include <fstream>
#include <cstring>
#include <unordered_set>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#include "evaluate.h"

void process_mem_usage(double& vm_usage, double& resident_set){
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}


void usage() {
    cerr << "./PRSim"<<endl
	 << "-d <dataset>"<<endl
	 << "-f <filelabel>"<<endl
	 << "-algo <algorithm>"<<endl
	 << "[-e <epsilon> (default 0.1)]"<<endl
 	 << "[-c <damping factor> (default 0.6)]"<<endl
	 << "[-qn <querynum> (default 50)]"<<endl
	 << "[-use_idx <whether use index> (default true,indexnum=sqrt(n))]"<<endl
	 << "[-check_dup <check duplicate edges> (default 0)]" << endl;
}


//Check parameters
int check_inc(int i, int max) {
    if (i == max) {
        usage();
        exit(1);
    }
    return i + 1;
}


bool cmp(const pair<int, double>& a, const pair<int, double>& b){
    return a.second > b.second;
}


void checkDuplicateEdge(string filename,string outFilename){
    int n;
    ifstream infile(filename);
    ofstream outfile(outFilename);
    infile >> n;
    outfile << n;
    cout<<"n="<<n<<endl;
    string line;
    unordered_set<string> edgeSets;
    int count = 0;
    while(getline(infile,line)){
        if(edgeSets.find(line) != edgeSets.end()){
            cout << "duplicate: " << line << endl;
        }
        else
	    outfile<<line<<"\n";
            edgeSets.insert(line);
        //if(count++ % 100000 == 0)
        //    cout << count << endl;
    }
    infile.close();
    outfile.close();
    cout << "Duplication checks done!" << endl;
}


int main(int argc, char **argv){
    int i = 1;
    char *endptr;
    string filename;
    string outputname = "-1";
    string filelabel;
    string algo;
    int querynum = -1;
    double eps = 0.1;
    double rmax = 0.0002;
    double c = 0.6;
    int k = 50;
    int indexnum=0;
    int isShowResult = 1;
    bool use_index = true;
    bool check_DuplicateEdge = false;
      
    if(argc < 4){
        usage();
        exit(1);
    }
    while (i < argc) {
        if (!strcmp(argv[i], "-d")) {
            i = check_inc(i, argc);
            filename = argv[i];
        } else if (!strcmp(argv[i], "-f")) {
            i = check_inc(i, argc);
            filelabel = argv[i];
        } else if (!strcmp(argv[i], "-algo")) {
            i = check_inc(i, argc);
            algo = argv[i];
        } else if (!strcmp(argv[i], "-e")) {
            i = check_inc(i, argc);
            eps = strtod(argv[i], &endptr);
            if ((eps == 0 || eps > 1) && endptr) {
                cerr << "Invalid eps argument" << endl;
                exit(1);
            }
	    rmax = 0.002 * eps;
        } else if (!strcmp(argv[i], "-c")) {
            i = check_inc(i, argc);
            c = strtod(argv[i], &endptr);
            if ((c == 0 || c > 1) && endptr) {
                cerr << "Invalid c argument" << endl;
                exit(1);
            }
        } else if (!strcmp(argv[i], "-qn")) {
            i = check_inc(i, argc);
            querynum = strtod(argv[i], &endptr);
            if ((querynum < 0) && endptr) {
                cerr << "Invalid querynum argument" << endl;
                exit(1);
            }
        } else if (!strcmp(argv[i], "-o")) {
            i = check_inc(i, argc);
            outputname = argv[i];
        } else if(!strcmp(argv[i], "-use_idx")) {
	    i = check_inc(i, argc);
	    int use_idx_flag = strtod(argv[i],&endptr);
	    if ((use_idx_flag!=0) && (use_idx_flag!=1)){
	      	cerr << "Invalid argument: -use_idx"<< endl;
		exit(1);
	    }
	    else if (use_idx_flag==0){
	    	use_index = false;
	    }
        } else if(!strcmp(argv[i], "-check_dup")){
	    i = check_inc(i, argc);
	    int check_dup_flag = strtod(argv[i],&endptr);
	    if ((check_dup_flag!=0) && (check_dup_flag!=1)){
	    	cerr << "Invalid argument: -check_dup"<< endl;
		exit(1);
	    }
	    else if (check_dup_flag==1){
	    	check_DuplicateEdge = true;
	    }
	} else {
            usage();
            exit(1);
        }
        i++;
    }

    if(filename==""|| filelabel=="" || algo==""){
	cout<<"Miss parameters"<<endl; 
   	usage();
	exit(1);
    }

    if(check_DuplicateEdge){
    	checkDuplicateEdge(filename,"dataset/"+filelabel+"_rm.txt");
	return 0;   
    }

    SimStruct sim = SimStruct(filename, filelabel, eps, rmax, c);
    if(use_index){
   	indexnum=5*sqrt(sim.vert);//index_num=j0 
    }
    
    if(querynum == -1 || querynum > sim.vert)
        querynum = min(50,sim.vert);
    

    if(algo == "SORT_INDEGREE_ASC"){
        sim.generateSortByIndegreeAsc("dataset/" + filelabel + "_sorted.txt");
    }


    //Output des-sorted pagerank results into the outfile
    if(algo == "GEN_PAGERANK"){
        sim.PageRankPowerMethod(50, "pagerank/" + filelabel + ".txt");
    }

    if(algo == "GEN_QUERY"){
        ofstream data_idx("query/" + filelabel + ".query");
        for(int i = 0; i < querynum; i++){
            data_idx << sim.R.generateRandom() % sim.vert << "\n";
        }
        data_idx.close();
    }

    //generate index,sort index and write index
    if(algo == "PREPROCESS"){
        cout<<"eps="<<eps<<endl;
	cout<<"rmax="<<100*rmax<<endl;
	cout<<"indexnum="<<indexnum<<endl;
	
	string pgname = "pagerank/" + filelabel + ".txt";
	ifstream input;
	input.open(pgname.data());
	if(!input.is_open()){
		cout<<"calculate pagerank"<<endl;
        	clock_t pg_start = clock();
		sim.PageRankPowerMethod(50,"pagerank/"+filelabel+".txt");
		clock_t pg_end = clock();
		cout << "calculating pagerank's time: "<< (pg_end - pg_start) / (double) CLOCKS_PER_SEC <<" s"<<endl;
		cout<<"====calculate pagerank done!===="<<endl;
		input.open(pgname.data());
	}
	assert(input.is_open());
	
	stringstream ss_idxname;
	ss_idxname<<"index/"<<filelabel<<"_"<<eps<<"_"<<indexnum<<".idx";
	ofstream data_idx(ss_idxname.str());

	double core_time=0.0;
        clock_t start = clock();
	for(int i = 0; i < indexnum; i++){
	    int nodeId;
            double tempSim;
            if(!(input >> nodeId >> tempSim)){
		break;
	    }
	    sim.backwardSearch(nodeId, data_idx, 100*rmax);
        }
	
        clock_t end = clock();

	for(auto &w: sim.idx_map){
            for(auto &level: sim.idx_map[w.first]){
                for(int i = 0; i < sim.idx_map[w.first][level.first].size(); i++){
                    data_idx << sim.idx_map[w.first][level.first][i].first << " " << w.first << " " << level.first << " " << sim.idx_map[w.first][level.first][i].second << endl;
                }
            }
        }

	cout << "Preprocessing time: "<< (end - start) / (double) CLOCKS_PER_SEC <<" s"<<endl;
	cout << "====preprocess done!===="<<endl;
	data_idx.close();
        input.close();
    }


    if(algo == "QUERY"){
	cout<<"eps="<<eps<<endl;
	cout<<"rmax="<<rmax<<endl;
	cout<<"indexnum="<<indexnum<<endl;
        
	double vm, rss;
        process_mem_usage(vm, rss);
        double lastRss = rss;

	stringstream ss_idxname;
	ss_idxname<<"index/"<<filelabel<<"_"<<eps<<"_"<<indexnum<<".idx";
        sim.loadIndex(ss_idxname.str());

	process_mem_usage(vm, rss);
        double indexMem = rss - lastRss;

        string pagerankname = "pagerank/" + filelabel + ".txt";
        ifstream fpagerank;
	fpagerank.open(pagerankname.data());
	assert(fpagerank.is_open());

        for(int i=0;i<indexnum;i++){
	    int tempNode;
            double tempSim;
            if(!(fpagerank >> tempNode >> tempSim)){
	    	break;
	    }
            sim.idx_node.push_back(tempNode);
        }
        
        string queryname = "query/" + filelabel + ".query";
	ifstream fquery;
	fquery.open(queryname.data());
	assert(fquery.is_open());

	cout<<"querynum="<<querynum<<endl;
        for(int i = 0; i < querynum; i++){
	    int nodeId;
            fquery >> nodeId;
            cout<<i<<": "<<nodeId<<endl;
            sim.singleSourceQuery(nodeId);
        }
    
	fquery.close();    
        cout << "index memory: "<< (indexMem/1024/1024)<<" GB" << endl;
	cout << "avg time: " << sim.avg_time / (double) querynum << endl;
	cout << "====single-source query done!===="<<endl;
    }
    
    return 0;

}

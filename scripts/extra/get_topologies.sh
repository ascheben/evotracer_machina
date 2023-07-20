python print_seeding_topology.py all_results_migration_optimal_split_22082022.txt | grep seeding| sed 's/ .*{/ /'| sed 's/}//'| cut -d' ' -f1,4,7,10,12,14,17 | sed 's/,//g' | tr ' ' '\t'

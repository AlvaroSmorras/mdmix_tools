#!/bin/bash
# @ job_name = crd
# @ initialdir = .
# @ wall_clock_limit = 40:00:00
# @ output = crd.q.o
# @ error =  crd.q.e
# @ total_tasks = 1
# @ tasks_per_node = 1
# @ cpus_per_task = 1

module load amber

for i in {1..100}; do
   if test -f md${i}.nc; then
      sed 's/XX/'$i'/g' tmp.traj > nc${i}crd.traj
      cpptraj -i nc${i}crd.traj
   fi
done

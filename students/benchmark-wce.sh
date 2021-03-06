#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
run_ce_solver()
{
	PROGRAMM_NAME="$1"
	LOG=$2
	CSV=$3
	maxSec=$4
	maxSecPerInstance=$5
	maxNotSolved=$6
	
	
	FILES=$(ls $7*.dimacs)

	rm -f time.txt

	overallTime=$(date +%s);
	now=$(date +%s);
	elapsed=`expr $now - $overallTime`;
	notSolved=0
	showedHint=0

	maxSecPerInstanceHard=$((maxSecPerInstance + 2));

	for f in $FILES
	do
		if [ $elapsed -le $maxSec -a $notSolved -le $maxNotSolved ]; then
			echo $f >> $LOG
			
			# start everything in a new process group such that we can kill everything if necessary
# 			(setsid /usr/bin/time -f "%e" -a -o time.txt timeout -k $maxSecPerInstance -s 9 $maxSecPerInstance $PROGRAMM_NAME< $f 1> prog_out.txt 2>&1) & PID=$!
			(setsid /usr/bin/time -f "%e" -a -o time.txt timeout --preserve-status -k $maxSecPerInstanceHard -s 2 $maxSecPerInstance $PROGRAMM_NAME< $f 1> prog_out.txt 2>&1) & PID=$!


			# kill processes when exiting this script
			trap "{ kill -$PID 2>/dev/null; }" TERM
			trap "{ kill -9 -$PID 2>/dev/null; }" EXIT

			wait $PID

			# just to be sure: if the process still is around brutally kill it
			kill -0 $PID 2>/dev/null || kill -9 -$PID 2>/dev/null;

			# get n
			n=$(grep "^[^#]" $f | head -1 | sed 's/#//' | sed 's/ .*//')

			# get m
			m=$(grep "^[^#]" $f | grep -v - | grep " " | wc -l)
			
			# get k
			#k=$(grep -ve "^#" prog_out.txt | wc -l)
			recursiveSteps=$(grep -e "#recursive steps:" prog_out.txt | sed -e 's/.*recursive steps: \([0-9]*\).*/\1/' )
			lastK=$(grep -e "last-k:" prog_out.txt | sed -e 's/.*last-k: \([0-9]*\).*/\1/' )
			cat prog_out.txt >> $LOG

			solveReason=$(grep -e "Finished:" prog_out.txt | sed -e 's/.*Finished: \([a-zA-Z]*\).*/\1/' )
			
			# get time
			time=$(cat time.txt);
			
			if [[ $time == "Command terminated by signal 9"* ]] || [[ $time == "Command exited with non-zero status"* ]]; then
				finished=0;
				(( notSolved += 1 ));
				time="";
			else
				finished=1;
			fi
			
			verify="";
			
			ourk=$(cat $f.solution); 
			if [ "$finished" -eq "1" ]; then
				solFile=$(basename $f .dimacs)
				# the best cost that we found
				msg=$($DIR/wce-verify.py $f prog_out.txt)
				if [ $? != 0 ]; then
					verify="verifier.py:  \t $msg";
				else
					k=$msg;
					if [ -n "$ourk" ] && [ "$ourk" -eq "$ourk" ] 2>/dev/null; then
						if [ "$ourk" -eq "$k" ]; then
							verify="correct";
						elif [ "$ourk" -gt "$k" ]; then
							verify=">>GOOD COST<<\tvs $ourk";
						else
							verify=">>BAD COST<< \t$(($k-$ourk))";
						fi
					else
						verify="our own cost is weird"
					fi
				fi
			fi
			
			fileNameLong=$(printf '%-40s' "$f");
			
			echo -e "$fileNameLong\t$n\t$m\t"$time"\t"$k"\t"$recursiveSteps"\t"$finished"\t"$verify"\t"$lastK"\t"$ourk"\t"$solveReason
			echo -e $f";"$n";"$m";"$time";"$k";"$recursiveSteps";"$finished";"$verify";"$lastK";"$ourk >> $CSV

			echo "" >> $LOG
			rm -f prog_out.txt
			rm -f time.txt

			now=$(date +%s);
			elapsed=`expr $now - $overallTime`;
		else
			if [ $showedHint -eq 0 ]; then
				if [ $notSolved -ge $maxNotSolved ]; then
					echo "$notSolved instances not solved. Script aborted."
				else
					echo "maximal time of $maxSec sec elapsed. Script aborted."
				fi
				showedHint=1;
			fi		
		fi
	done
}


j="julia"

if [ $PWD == "/home/users/h/hjalmar.schulz/Desktop/AlgEng" ]; then
	echo "Using uni julia"
	j="/home/users/h/hjalmar.schulz/julia"
fi


PROGRAMM_NAME=$j" ./julia/main.jl"  		# insert your program here

today=$(date +%Y-%m-%d-%H-%M-%S)

LOG="log.txt"								# specify the name of the log file
maxSec=432000								# overall allowed time for the whole script
maxSecPerInstance=300							# allowed time (in seconds) for one instance
maxNotSolved=20000								# no of instances the program is allowed to fail to solve. If reached, then the script is aborted

CSV=$1".csv"
if [ -z "$1" ]; then
	CSV="results-"$today".csv"
fi
echo "file;n;m;time;solsize;recsteps;finished;verified;lastK;opt" > $CSV
## now loop through data set directories
for data in $(find $DIR -mindepth 1 -maxdepth 1 -type d); do
	FILENAME=$(basename $data)
	echo "run $data instances $PROGRAMM_NAME (Tab-separated columns: File, n, m, Time in seconds, solution size, recursive steps, finished, solution size verified)"
	run_ce_solver "$PROGRAMM_NAME" $LOG $CSV $maxSec $maxSecPerInstance $maxNotSolved $data/
done

echo ""

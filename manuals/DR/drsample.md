# Configurable parameters STDR

* Number of nodes to use per replica
* Number of backups to make
* Number of steps to run per backup
* Check how much time left, if less than T amount of time, backup

- How is DRPE computed if a single replica crashes?
	- take the stationary position of the crashed replica and use that to compute DRPE(?)
	
		

# DR Feature requests 
* all sending of non-numerical data as additional data for server


# June 27 2010

1) Expected 131 x 100 moves in totoal but got 28,000 moves 
	- is this because DR keeps going if the cpu doesn’t shut down?
2)  How to interpret the analyzeforcedb graphs?
3) Fix failures - this is because too many backups causing the thing to fail
4)  Compare replica sampling density with Sarah’s 
5) Do not allow multiple replica @ same nominal position
	- makes sense now that I think about it
	- It is really a setting for temperature range and cancellation energies table for STDR run 
6) If you’re restarting a simulation using STDR, can encode the temperature information in gromacs specific files
    - what is the snapshot file?

7) Document error messages from the DR package!!!!

Error: replica temperature must be unique and go in descending order
 - cannot have repeat nominal positions

Error: the force database header has inapropriate parameters
- cannot change the number of additional data once database is started and if you want to continue to write to the forcedb

8) Allow sending of data that is non-numerical to the server and store this
	- allow polymorphism in data type - how to do this properly in C?
	

# July 2nd 2010 

* Stuff I learned about DR package from playing around with it
	* Tricking STDR running N replicas at hard coded temps and writing out a snapshot
	* This doesn’t work because the server never knows what temperature that a client actually ran at.  It assumes that the client ran at the temperature or whatever setting it passed the client. (Assumes that client is well behaved)
	* Reassign temps & restarting uniformly 
		* I don’t want this because this is what I think is going on:
			- Simulating a small segment at a different temperature means reequilibrating system at a new temperature.  And if the system from that small segment is not representative of the system at temperature T, a subsequent move from that temperature to another temperature is equivalent of the improbably move between two temperatures with non-overlapping potential energy distributions.  --- I don’t know enough of this physics to say -- but I think it is unphysical and will introduce artifacts, though over time the algorithm (STDR) should correct this (just think of how the simulation was started ... each structure was simulated at 800K, then simulated at their starting T for 9 ps and requires a long simulation time to “recover”)

## Notes on current implementation:
	
DR_server.cpp:
In main(), a single object of  class ‘read_input_script_file_class’ calls input_script->read_input_script_file(..,..)

This function appears to parse the ‘.script’ file and sets the initial states of all replicas using the JOB rows as specified in ‘.script’.  Also parses the KEYWORDS (see documentation).  Effectively initialize the

# DR caveats:
* Issue: When restarting from a snapshot.  The server quits without telling why.
Explanation:  This was caused by the .script file having the same number of completed steps as written in the snapshot file.  Server sees this and thinks simulation should end and quits.
* Fix: Should output a warning / info type message to user so user knows to correct the script file
* Simple test case to reproduce: 
	* Specify 2 replicas with a single move
	* Save a snapshot
	* Restart system without changing the script file

# Test result:
* log file will contain “####### Session Ended ####” as the last line with no other messages indicating why server quits.


# DR improvement ideas:
* Refactor DR so that the code responsible for VRE/ST/DR/MC is separate from code dealing with server/client communication and utility code for initializing server, replica, etc. This way, DR_server.cpp will be short and therefore easier to read, understand and modify

* Documentation: Note that runs starting with differing additional data should make sure the older ‘forcedatabase’ file is not present so a new one is written otherwise, server will quit with ‘header error’

* Clarify how a ‘.snapshot’ file is utilized in restarts.  Clarify what exactly is encoded in the ‘.snapshot’ file

# More DR problems that I've noticed:

- additional data ... if still in the file is not a float .. it doesn’t check the integrity of the data ... and just puts in garbage ...
- should be allowed to characters/ strings as well as numbers ...
- not designed to send a lot of data 

- rewrite code to take multi column file and write it to the database

ie 
analysis writes 
file 1: with A B C
file 2 with D E F
should just be written as

A B C D E F 


- when does the server output steps ?
   -- sometimes a run will write to database but not to the output directories ... why?

- check the variable CPTtime at the top of the DR_client_wrapper

- error from analyzeforcedatabase

calcsample() is zero ... 
data text file is still outputted


Annoyances

- sending add data when no add data keyword is added to teh script ... the data is not written to db ... nothing is written to db


# TODO (C. Neale)
* Important:
	* The snapshot should really keep a record of Nsamesystem_uncoupled and vre secondary and primary system sizes
	* Go through what Tom has suggested in read_input_script_file.h after reading the script and perhaps allow some disallowed options
	* it is unclear why the tester does not exit when I give the command via '-r' to use fewer replicas
	* Go through the source and see where things like Temperature have been left to an expective else statement especially for move type as as there is also the NoMoves statement and it will reduce future errors if I introduce new types of moves or rx coords

For Later:
* analyze_force_database gives me a segmentation fault for NoMoves:
                * Reading data for force plots
                * Number of discrete w positions: 1
                * max time is: 1998
                * Segmentation fault (core dumped)

* Would be nice to allow automation of multiple rounds of cancellation.  However, I am moving towards using vRE, so cancellation is not important to CN any more.

* fix the nligands routines
   PAUSED, because:
     look at this from read_input_script_file.h:
       else if(strcasecmp(param,"LIGAND2")==0) column[i]=ColumnW2;
       else if(strcasecmp(param,"FUNNEL")==0) column[i]=ColumnW2;
     So I need to talk to Rowan and Tom about what they really mean. The script file indicates that they have
     entirely different uses, but now it does not appear to be actually coded that way.

* test that different force constants for different restraints are set up correctly for umbrella while using an exchange method other than non-monteCarlo

* make the circular coordinate read-in a little easier to use

* The sampling density is displayed as 1/2 when nligands = 2.
     - What other things might this influence?

* add ticks at the nominal positions of plot

* For the tester: print out the cancellation values when loading a snapshot

* Add WHAM solution to analyze_force_database

* Add ability to calculate the cancellation values via WHAM instead of simple force integration

* Add ability for analyze_force_database to calculate the MSD or some other measure of motion along the rx. coord

	* move the standard definitions (e.g. Boltzmann's constant) to a single header = defs.h
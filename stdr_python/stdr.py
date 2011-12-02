
Job
- abstraction for a job running on a cluster node
- contains a replica
- and information to dispatch a replica

Job db



Replica
    - represents a markov state
    - with temperature, sequence number
    - knows the next move?

Replica db - I think keeping a separate db for replica is redundant.  Job is a replica but with additional node information.
We want to ensure that if Job data is inserted, replica is inserted
There is a one to one association with Job - Replica, replica is dependent on Job or Job is dependent on Replica?


every monte carlo step produces
sequence number (monte carlo step number)
lambda (or temperature in ST's case)
    - beta
    - replica position as an integer i = 0 and N-1, where N is the total number of replicas


server:
- qsub N number of clients (requesters)
- initial step
    - when a client starts up 
    - makes a GET request to the server for a Job object
    - the server determines the simulation parameters to be sent back and responds serially for each replica  
    - after some time, the client sents back the Job object with a PUT with the data to be recorded in the database
        - reports status of the simulation: STOPPED, ERROR, TIMEDOUT, SUCCESS
    - server serially accesses db and records job information 

- store list of replica positions (or current temperatures)
- implements probability calculation (transition probability for a replica) and calculation of the DRPE
- access the database
    - writes to the replica parameter database the last simulation that was finished
    - determine parameters for the next simulation (responds with a new Job request object)




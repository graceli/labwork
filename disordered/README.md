# Notes
* This is a note

# TODO
* Refactor out the code relating to job submission and state checking to a separate module
	* Keep a list of files that have been processed.  All files are unprocessed initially
	* have one script that reports on the state of the analysis
		* show failed files and their error messages
	
	* Option to allow resubmission of failed analysis jobs
			* resumes and reanalyze failed analyses
			* report on the number of failed trajectories, eg:
					- "analyzing 120 trajectories ... "

* Define an interface to this package
    * Allow users to extend functionality by subclassing
        * Add customized analysis
        * Add new filetype support
        
* Separate preprocessing and ETL code from the analysis code

## List of Modules
* submission
* etl
* analysis
    * convergence


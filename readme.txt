Hi Dr. Ridley!

First in the defualts.json, add this on line 59 
"Logfile" : {"name" : "UA/output/log.txt"},

Then, inside outputs in the defaults.json add 
	"species" : ["O2", "O2+"]	

Now try to run the program. There should be an error with the logfile in run.test UA/output/log.txt such that nothing shows up.
I am not sure what the error is in logfile.cpp that is causing this. I believe there is a for loop that is not exiting, but I am
not sure because if I add cout statements they show up. Let me know what you think.

-Rishi
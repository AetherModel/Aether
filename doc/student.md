# Student Walk Through of Aether


## Start with the Basics

0. [Using the terminal - Unix commands](https://ubuntu.com/tutorials/command-line-for-beginners#1-overview)

1. [A good tutorial on cmake](https://www.internalpointers.com/post/modern-cmake-beginner-introduction)

2. [Information on the apt package manager](https://ubuntu.com/server/docs/package-management)

3. [Information on the port package manager](https://brianreiter.org/2020/06/25/why-and-how-to-set-up-macports-package-manager-for-macos/), [and a youtube video on this topic](https://www.youtube.com/watch?v=fYY7mArCryI)

4. [Tutorial for git and GitHub](https://www.freecodecamp.org/news/git-and-github-for-beginners/)


## Some Editors

The Classics:

1. [Emacs](http://www.jesshamrick.com/2012/09/10/absolute-beginners-guide-to-emacs/)

emacs can be installed using a package manager such as apt or port

2. [vim](https://linuxconfig.org/vim-tutorial)

21st Century Editors:

1. There are always popular Integrated Development Environments (IDEs) that can be used to 
    edit, develop, and test code.  If you use one and develop on Aether, please feel free to update
    this section.


## Install Aether for Development

Make sure you have all of the dependencies installed (see README.md).

1. Create a new fork of the Aether repository on github.  To do this:

1a. Go to the Aether repository and click on "Fork" then "+ Create a new fork".

1b. Unclick the "Copy the main branch only" checkbox

1c. Click "Create fork"

1d. github will take you to this repository automatically.  It is now
your "version" of Aether, and all development will be done here.  I am going to use the username "YourRepo" instead of your actual repository below.


```bash
git clone https://github.com/YourRepo/Aether
cd Aether
git checkout develop
```

To compile Aether, you need to make sure you have all of the dependencies
installed (see the main README.md in the Aether directory), and then run
these commands:
```bash
mkdir build
cd build
cmake ..
make
```

Once you have compiled you can run Aether with the standard run directory
like this:
```bash
cd ..
cp -R share/run ./run.test
cd run.test
./aether
```

You should see something like:
```bash
run.test% ./aether
> Need to NOT adjust F10.7, but that isn't included yet!!!
> Writing file : 3DALL_20110320_000000
> Writing file : 3DBFI_20110320_000000
> Wall Time : 4s (left : 1111h); Current Time : 2011 3 20 0 0 0 0 
> Wall Time : 4s (left : 23m); Current Time : 2011 3 20 0 0 10 0 
> Wall Time : 5s (left : 14m); Current Time : 2011 3 20 0 0 20 0 
> Wall Time : 5s (left : 9m); Current Time : 2011 3 20 0 0 30 0 
> Wall Time : 5s (left : 7m); Current Time : 2011 3 20 0 0 40 0 
> Wall Time : 6s (left : 7m); Current Time : 2011 3 20 0 0 50 0 
> Wall Time : 6s (left : 5m); Current Time : 2011 3 20 0 1 0 0 
> Wall Time : 6s (left : 5m); Current Time : 2011 3 20 0 1 10 0 
> Wall Time : 7s (left : 5m); Current Time : 2011 3 20 0 1 20 0 
> Wall Time : 7s (left : 4m); Current Time : 2011 3 20 0 1 30 0 
etc.
```

Now, edit the aether.json file and add the following after the "Debug"
section, where "Your Name" should be your name:
```bash
    "Student" : {
	"name" : "Your Name",
	"is" : true },
```

Then, run aether again.  Your output should now look something like this:
```bash
run.test% ./aether
> Hello Aaron - welcome to Aether!
> Need to NOT adjust F10.7, but that isn't included yet!!!
> (4) What function is this Aaron?
> (3) What function is this Aaron?
> Writing file : 3DALL_20110320_000000
> Writing file : 3DBFI_20110320_000000
> Wall Time : 4s (left : 1111h); Current Time : 2011 3 20 0 0 0 0 
> (1) What function is this Aaron?
> (2) What function is this Aaron?
> (3) What function is this Aaron?
> (1) What function is this Aaron?
> (3) What function is this Aaron?
> Wall Time : 4s (left : 23m); Current Time : 2011 3 20 0 0 10 0 
> (1) What function is this Aaron?
> (3) What function is this Aaron?
> (1) What function is this Aaron?
> (3) What function is this Aaron?
> Wall Time : 5s (left : 14m); Current Time : 2011 3 20 0 0 20 0 
```

There are 5 messages in this output (4 "what is this function"
messages and a "welcome" message). Now that you have this working, you
can try to do some development on the code.  The first step in doing
development is making a new branch.  To do that, you need to run a
command like this:

```bash
git checkout -b my_new_branch
```

where my_new_branch should be named something that makes sense.

Once you have done that, you can start editing files. The goal of this
part is to get you familiar with editing files and finding things in
the code. Find all of the messages and add some sort of response to
all five places in the code where these messages occur.  For example:

```bash
run.test% ./aether
> Hello Aaron - welcome to Aether!  Thanks! 
> Need to NOT adjust F10.7, but that isn't included yet!!!
> (4) What function is this Aaron? Found it: xyz
> (3) What function is this Aaron? Found it: abc
> Writing file : 3DALL_20110320_000000
> Writing file : 3DBFI_20110320_000000
> Wall Time : 4s (left : 1111h); Current Time : 2011 3 20 0 0 0 0 
> (1) What function is this Aaron? Found it: abc
> (2) What function is this Aaron? Found it: def
> (3) What function is this Aaron? Found it: ghi
> (1) What function is this Aaron? Found it: abc
> (3) What function is this Aaron? Found it: ghi
> Wall Time : 4s (left : 23m); Current Time : 2011 3 20 0 0 10 0 
> (1) What function is this Aaron? Found it: abc
> (3) What function is this Aaron? Found it: ghi
> (1) What function is this Aaron? Found it: abc
> (3) What function is this Aaron? Found it: ghi
> Wall Time : 5s (left : 14m); Current Time : 2011 3 20 0 0 20 0 
```

Once you have found all of them, then altered the code to note that
you have found them, you can then commit the code and then push the
code to your repository.  When that is all done, send a message to
your advisor, letting them know that you have completed this task.

While you are awaiting a response from your advisor, you can check out
the Issues part of the Aether github site.  That will let you know all
of the outstanding issues that are being worked on. Maybe you can
contribute to some of these!


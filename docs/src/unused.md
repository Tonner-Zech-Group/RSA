## Using Multithreading with VSCode
Running RSA simulations over VSCode in parallel is possible by changing the "julia.NumThreads" setting. Simply search for "num thread" in the settings searchbar and change to your needs. 
To test the new settings use the following command in your notebook:
```
Threads.nthreads()
```
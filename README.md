# Limit Scan

To create datacards for the limit scan wrt BDT cuts, use:

~~~
./runToScan.py
~~~

Before running using `--pdf_type=flat`, run using `--pdf_type=unfixed_exp` first.

Once the datacards, are made, run:

~~~
./readLimit.py
~~~

for submitting the bdt scan runs to condor and for plotting.

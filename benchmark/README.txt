Files used:

code_generator/ ( Local Computer )
  Algorihtm.rb : to generate the code that is outputed to CallTemplate.f90 and MagicFilter.f90. The generated files should be put in the bigdft-1.6.0/src/convolutions-c/ folder whereafter conv_check is compiled with "sh run.sh -argument" 
  CallTemplate.f90 : Contains the code for starting the starting the counters, call the convolution function, stop and output the result of the counter.
  MacigFilter.f90 : Contains the generated convolution functions.

bigdft-1.6.0/src/convolutions-c/ ( Tibidabo )
  batch_run.sh : called with the qsub command to run on a free node, calls run.sh.
  run.sh : runs the benchmark if no arguments are given. if any argument is provided it will compile each of the conv_check versions, one for every counter. each compiled binary is put in bin/.

  bin/ : contains the different versions of conv_check, one for every counter type.
  data/ : contains the resulting output of the benchmark.
  tests/ : contains the specific papi_counters.f90 files for each counter test

  conv_check.f90 : modified to include the code
                  INCLUDE 'papi_counters.f90'
                  INCLUDE 'CallTemplate.f90'
                  INCLUDE 'MagicFilter.f90'

  Makefile/.in/.am : modified so that MagicFilter_sse.c is not used.
  extract_data.py : extracts the values of each test and puts them in the extracted_data directory
  fort.1 : contains the argument to be passed to the conv_check

analysis/ ( Local Computer )
  get_data_and_analyse.sh : downloads the extracted data from the server and runs present_results.R
  present_results.R : reads the extracted data and draws graphs
  extracted_data/ : where the downloaded data will be saved


Workflow:

  Local Computer:
    ruby Algorihtm.rb # Outputs the generated files

  Tibidabo:
    Put CallTemplate.f90 and MagicFilter.f90 in bigdft-1.6.0/src/convolutions-c/
    sh run.sh asdf # Compiles the 6 versions of the conv_check
    qsub batch_run.sh # Runs the tests on a free node. Note: The progress may be seen by tailing the stdout textfile that is created in the home directory folder, name of the file is convolution-test.o1234. May be edited in the batch_run.sh file.

    python extract_data.py # Extracts the data from the resulting tests

  Local Computer:
    sh get_data_and_analyse.sh # Rplots.pdf is created with the graphs

For large data sets the tests might take several minutes to complete. If the time set in batch_run.sh is not enough ( 25 minutes ) it may be increased.

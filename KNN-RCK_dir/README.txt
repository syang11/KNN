KNN-RCK extends RCK to accept a KNN predicted PWM as input. It then converts the input PWM to the k-mer sequence model by scoring each k-mer using the PWM, and applies this model in the binding preference prediction. Please refer to the manuscript and its supplementary note for details.

1. Here, the KNN-RCK file is an executable that extends RCK. The training process runs similiarly as RCK, i.e. most command line options passing to the programs are the same except for:
 -w    <motifwidth range> (eg. 7-7)
          controls the size range of the motif; since the PWM width is fixed in KNN-RCK, the two numbers need to be same: eg. when running the program with -w 7-7, the algorithm uses a input PWM with width equal to 7 and searches for motifs starting from width 7 until width 7.
 -f    <PWM file name> (eg. vts1.pwm)
          takes a PWM file as input; here the PWM file has to be in the format of a matrix with rows standing for {A,C,G,U} and columns standing for position 1 to 7.

Example Run:
Training (use the same example as in the RCK package [1]. Note here we have "-w 7-7" and "-f vts1.pwm"):
./KNN-RCK -b 200 -w 7-7 -a ACGU -e PHIME -f vts1.pwm -s 3 -c VTS1_training_sequences.txt -h VTS1_training_annotations.txt -d VTS1_test_sequences.txt -n VTS1_test_annotations.txt -o VTS1_demo -m ./outputs/

Predicting (same as RCK):
./KNN-RCK -w 4-4 -a ACGU -e PHIME -d VTS1_test_sequences.txt -n VTS1_test_annotations.txt -l VTS1_demo -m ./outputs/

2. To rebuild this executable, just replace the rnacontext.cpp file in RCK package [1] with the rnacontext.cpp file provided here. Then follow the same "make" process for RCK. 

[1] The original RCK program can be downloaded from http://cb.csail.mit.edu/cb/rck/

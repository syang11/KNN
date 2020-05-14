KNNpwm.pl requires 4 argument: <homolog PWM dir with PWMs having equal length> <protein similarity matrix file derived from ClustalW> <optK of KNN> <output dir for predicted PWMs>
  
* To predict PWM based on homolog neighbors' PWMs:
```
perl KNNpwm.pl homologs_PWM_dir protein_PID_file optK output_dir
```

* The `protein_PID_file` above is obtained from:
```
perl makePromat.pl input_dir fasta_file_name output_dir protein_PID_file
```
The `fasta_file_name` here can be sample.fa for a simple example.

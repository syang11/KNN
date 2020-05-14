KNNpwm.pl requires 4 argument: 
1. homolog PWM dir with PWMs having equal length
2. protein similarity matrix file derived from ClustalW
3. optK of KNN
4. output dir for predicted PWMs>
  
* To predict PWM based on homolog neighbors' PWMs:
```
perl KNNpwm.pl homologs_PWM_dir protein_PID_file optK output_dir
```

* The `protein_PID_file` above is obtained from:
```
perl makePromat.pl input_dir fasta_file_name output_dir protein_PID_file
```
The `fasta_file_name` here can be sample.fa for a simple example.

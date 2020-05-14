
* To predict PWM based on homolog neighbors' PWMs:
```
perl KNNpwm.pl homologs_PWM_dir proteinPIDfile optK output_dir
```

* The `proteinPIDfile` above is obtained from:
```
perl makePromat.pl input_dir fasta_file_name output_dir matrix_file_name
```

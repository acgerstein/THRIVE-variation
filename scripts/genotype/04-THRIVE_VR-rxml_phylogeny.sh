#custom example cripts for YST6

raxml-ng --all --msa 230623_YST6_fasta.min1.fasta  
--model GTR+G --prefix 
YST6 --threads 15 --seed 3 --bs-metric fbp,tbe

raxml-ng --all --msa TVY10.min1.fasta --model GTR+G --prefix TVY10 
--threads 15 --seed 15 --bs-metric fbp,tbe

raxml-ng --all --msa TVY4.min1.fasta --model GTR+G --prefix TVY4 
--threads 15 --seed 15 --bs-metric fbp,tbe

raxml-ng --all --msa TVY7.min1.fasta --model GTR+G --prefix YST7 
--threads 15 --seed 15 --bs-metric fbp,tbe

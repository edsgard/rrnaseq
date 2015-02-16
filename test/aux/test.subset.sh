

#Subset rows and cols
cd test
subset -i ./rqc/data/rpkm.rds -o ./rqc/data/rpkm.test_subset.rds -r ./subset/genes.txt -c ./subset/cells.txt
subset -i ./rqc/data/rpkm.rds -o ./rqc/data/rpkm.test_subset.x.rds -r ./subset/genes.txt -c ./subset/cells.txt -x
subset -i ./rqc/data/rpkm.rds -o ./rqc/data/rpkm.test_subset.e.rds -r ./subset/genes.txt -c ./subset/cells.txt -e
subset -i ./rqc/data/rpkm.rds -o ./rqc/data/rpkm.test_subset.xe.rds -r ./subset/genes.txt -c ./subset/cells.txt -x -e

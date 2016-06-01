

#Subset rows and cols
cd test

#subset rows
subset -i ./rqc/data/rpkm.rds -o ./rqc/data/rpkm.test_subset.rds -r ./subset/genes.txt

#subset rows and cols
subset -i ./rqc/data/rpkm.rds -o ./rqc/data/rpkm.test_subset.rds -r ./subset/genes.txt -c ./subset/cells.txt

#filter rows, subset cols
subset -i ./rqc/data/rpkm.rds -o ./rqc/data/rpkm.test_subset.x.rds -r ./subset/genes.txt -c ./subset/cells.txt -x

#subset rows, filter cols
subset -i ./rqc/data/rpkm.rds -o ./rqc/data/rpkm.test_subset.e.rds -r ./subset/genes.txt -c ./subset/cells.txt -e

#filter rows and cols
subset -i ./rqc/data/rpkm.rds -o ./rqc/data/rpkm.test_subset.xe.rds -r ./subset/genes.txt -c ./subset/cells.txt -x -e

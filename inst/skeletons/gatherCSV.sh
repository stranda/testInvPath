#!/bin/bash

rm -fr tmp.csv header.csv
#mv *.csv csv/
mv results/*.csv archivedCSV/
echo "moved csvs"
cat archivedCSV/*.csv >tmp.csv
head -n1 tmp.csv > header.csv
echo "made tmp"
cat tmp.csv |sed '/introNe/d' > tmp2.csv
cat header.csv tmp2.csv > gathered/`date|tr ' ' '_'|tr ':' '_'`.csv

#echo "made gathered object, tarrin up old csvs"

cd archivedCSV
#tar czf `date|tr ' ' '_'|tr ':' '_'`.tar.gz *.csv
rm *.csv
cd ..

echo "cleaning up"
rm tmp.csv header.csv tmp2.csv
echo "done"


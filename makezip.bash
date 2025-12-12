rm -r zip
mkdir zip
cd zip
cp -r ../include .
cp -r ../src .
cd ..
tar -czf ./md.tar.gz ./zip

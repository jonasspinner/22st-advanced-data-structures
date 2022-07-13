
for name in bv_10k\
            bv_100k\
            bv_1M\
            tree_d9_c4\
            tree_d6_c5-10\
            tree_d8_c5;
do
  echo "example_${name}.txt"
  curl --request GET -sL \
       --url "https://algo2.iti.kit.edu/download/example_${name}.txt"\
       --output "example_${name}.txt"
done

for name in bv_10k\
            bv_100k\
            bv_1M;
do
  echo "example_${name}_output.txt"
  curl --request GET -sL \
       --url "https://algo2.iti.kit.edu/download/example_${name}_output.txt"\
       --output "example_${name}_output.txt"
done
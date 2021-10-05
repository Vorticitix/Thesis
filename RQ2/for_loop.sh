arr1=($(seq 0 4 8))
arr2=($(seq 4 4 12))
arr3=($(seq 30 4 38))
arr4=($(seq 34 4 42))
echo "${arr2[@]}"
for i in $(seq 1 "${#arr1[@]}")
do 
	for j in $(seq 1 "${#arr3[@]}")
	do 
		echo "$i"
		echo "$j"
		echo "${arr1[i-1]}"
		echo "${arr2[i-1]-10}"
		echo "${arr3[j-1]}"
		echo "${arr4[j-1]}"
	done
done 
# for i in "${arr[@]}"

# do
# 	for j in "{i[@]}"
# 	do
# 		echo "$j"
# 	done
# done
# # echo ${arr1[@]}

# array1=(name1 name2)
# name1=(one two)
# name2=(red blue)


# for name in "${arr[@]}"
# do
#   typeset -n nameref="$name"
#   for value in "${nameref[@]}"
#   do
#     a="${value[0]}"
#     echo "$a"
#   done
# done

# for name in "${arr[@]}"
# do
# 	typeset -n nameref="$name"
# 	a="${nameref[1]}"
# 	echo "$a"
# done
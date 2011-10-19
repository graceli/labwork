#!/usr/bin/awk

{
	if($1 ~ /G/){sum+=$1};
	if($1 ~ /M/){sumM+=$1};
	if($1 ~ /k/){sumK+=$1}
}
END{
	print sum,sumM,sumK
}

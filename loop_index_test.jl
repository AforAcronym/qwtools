# loop_index_test
export test_idx 

function test_idx(rows::Int64, cols::Int64, depz::Int64)
	# cols  = 700;
	# rows  = 700;
	# depz = 700;

	print("rand(rows, cols, depz):  ")
	@time mx  = rand(rows, cols, depz);
	print("zeros(rows, cols, depz):  ")
	@time mx1 = zeros(rows, cols, depz);

	print("DRC: ")
	@time for d = 1:depz, r = 1:rows, c = 1:cols 
		mx1[r, c, d] = 1+mx[r, c, d] 
	end 

	print("DCR: ") # 2nd place
	@time for d = 1:depz, c = 1:cols, r = 1:rows 
		mx1[r, c, d] = 1+mx[r, c, d] 
	end 

	print("RDC: ")
	@time for r = 1:rows, d = 1:depz, c = 1:cols 
		mx1[r, c, d] = 1+mx[r, c, d] 
	end 

	print("CDR: ")
	@time for c = 1:cols, d = 1:depz, r = 1:rows 
		mx1[r, c, d] = 1+mx[r, c, d] 
	end 

	print("RCD: ")
	@time for r = 1:rows, c = 1:cols, d = 1:depz 
		mx1[r, c, d] = 1+mx[r, c, d] 
	end 

	print("CRD: ")
	@time for c = 1:cols, r = 1:rows, d = 1:depz 
		mx1[r, c, d] = 1+mx[r, c, d] 
	end 

	print("Length:  ") # Winner
	@time for i = 1:length(mx)
		mx1[i] = 1+mx[i] 
	end 

	print("Length + ind2sub: ")
	@time for i = 1:length(mx)
		r, c, d = ind2sub((rows,cols,depz), i)
		mx1[r, c, d] = 1+mx[r, c, d] 
	end 
end
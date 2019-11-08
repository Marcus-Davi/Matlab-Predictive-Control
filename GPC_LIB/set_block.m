% set_block(M,i,j,blksize,setblock)
% seta o bloco [i,j] de tamanho blksize da matriz M igual a setblock
function m = set_block(M,i,j,blksize,setblock)
%blksize -> [i j]
m = M;
%blksize -> [i j]
blk_row = blksize(1);
blk_col = blksize(2);

blki = (i-1)*blk_row+1;
blkj = (j-1)*blk_col+1;

m(blki:blki+blk_row-1,blkj:blkj+blk_col-1) = setblock;

end
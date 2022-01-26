## set up test matrix type 1
## m = (1 - tG)
## blockSize 
## numBlocks

l = 0
blockSize = 2*(l+1)^2
numBlocks = 3

n = blockSize * numBlocks

G0 = hilb(n);

t1 = hilb(blockSize)

t = t1;

for i = 2:numBlocks
    t = blkdiag(t,t1);
endfor
    
m = eye(n) - t*G0

t1 = postpad(t1,n,0,1)

t1 \ m


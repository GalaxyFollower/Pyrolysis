function aint=integratePolynomials(a)
[np, q]=size(a);
aint=[a./repmat((q:-1:1),[np 1]) zeros([np 1])];
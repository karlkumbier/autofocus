function templateImage = generateTemplate(a,b,c)

if nargin ==0, 
	a = 64; b = 32; c = 3;
end;
assert(a>=b);
templateImage = zeros(2 * b, 2 * a,c);
for j = 1 : a
	len = floor(b * sqrt(1-((j - 0.5)/a)^2) + 0.5);
	if (len > 0)
		templateImage(b + (1-len : len), a + [1-j;j],:) = 1;
	end
end; 

end

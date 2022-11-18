a = Channel.from('a','b','c')
b = Channel.from(1,2,3)
c = Channel.from('p','q')

c.concat( b, a ).view()

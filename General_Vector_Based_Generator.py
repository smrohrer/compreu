import string



ax=0
ay=0
az= 0

bx= 0
by=1.4
bz=0

v1x=2.425

v1y=0
v1z=0



v2x=1.226
v2y=2.12349
v2z=0


lax=1.2197
lay=-0.68727
laz=0
## la(z,y,z) are the base coordinates for the bottom zig-zag edge
lbz=0
## lbz is the z-coordinate for zig-zag edge of the top edge

name="Graphene"

m=7
n=7



print string.rjust("C",0), string.rjust(`lax`, 12), string.rjust(`lay`, 12), string.rjust(`laz`, 12)


oddvar=2

for i in range(m):
    for ii in range(n):
        if ii in range (0,2):
            t1x=ax+(i*v1x)+(ii*v2x)
            t1y=ay+(i*v1y)+(ii*v2y)
            t1z=az+(i*v1z)+(ii*v2z)
            
            t2x=bx+(i*v1x)+(ii*v2x)
            t2y=by+(i*v1y)+(ii*v2y)
            t2z=bz+(i*v1z)+(ii*v2z)

            rt1x=(round(t1x, 6))
            rt1y=(round(t1y, 6))
            rt1z=(round(t1z, 6))
        
            rt2x=(round(t2x, 6))
            rt2y=(round(t2y, 6))
            rt2z=(round(t2z, 6))

            print string.rjust("C",0), string.rjust(`rt1x`, 12), string.rjust(`rt1y`, 12), string.rjust(`rt1z`, 12)
            print string.rjust("C",0), string.rjust(`rt2x`, 12), string.rjust(`rt2y`, 12), string.rjust(`rt2z`, 12)

    
##These unit cells require no leftward movement corrections, & thus are run normally

                    
        elif ii in range (2,4):
            t1x=ax+(i*v1x-v1x)+(ii*v2x)
            t1y=ay+(i*v1y)+(ii*v2y)
            t1z=az+(i*v1z)+(ii*v2z)
            
            t2x=bx+(i*v1x-v1x)+(ii*v2x)
            t2y=by+(i*v1y)+(ii*v2y)
            t2z=bz+(i*v1z)+(ii*v2z)
            
            
            rt1x=(round(t1x, 6))
            rt1y=(round(t1y, 6))
            rt1z=(round(t1z, 6))
        
            rt2x=(round(t2x, 6))
            rt2y=(round(t2y, 6))
            rt2z=(round(t2z, 6))

            print string.rjust("C",0), string.rjust(`rt1x`, 12), string.rjust(`rt1y`, 12), string.rjust(`rt1z`, 12)
            print string.rjust("C",0), string.rjust(`rt2x`, 12), string.rjust(`rt2y`, 12), string.rjust(`rt2z`, 12)

            

    
##The 2nd & 3rd rows require only one leftward unit cell movement 


        elif (ii>3 and ii % 2 == 1):
            po=ii/2


            t1x=ax+(i*v1x)-(po*v1x)+(ii*v2x)
            t1y=ay+(i*v1y)+(ii*v2y)
            t1z=az+(i*v1z)+(ii*v2z)
                        

                        
            t2x=bx+(i*v1x)-(po*v1x)+(ii*v2x)
            t2y=by+(i*v1y)+(ii*v2y)
            t2z=bz+(i*v1z)+(ii*v2z)
                                       
                    
            rt1x=(round(t1x, 6))
            rt1y=(round(t1y, 6))
            rt1z=(round(t1z, 6))
        
            rt2x=(round(t2x, 6))
            rt2y=(round(t2y, 6))
            rt2z=(round(t2z, 6))

            print string.rjust("C",0), string.rjust(`rt1x`, 12), string.rjust(`rt1y`, 12), string.rjust(`rt1z`, 12)
            print string.rjust("C",0), string.rjust(`rt2x`, 12), string.rjust(`rt2y`, 12), string.rjust(`rt2z`, 12)

##^^^This will push all even rows above row 3 to the left to align in armchair with the previous rows since vectors t1 and t2 
##continually move the 2 carbon unite cell to the right, as they iterate through the loop generating the sheet




        elif (ii>3 and ii % 2 == 0):
	    pe=ii/2
            
            t1x=ax+(i*v1x)-(pe*v1x)+(ii*v2x)
            t1y=ay+(i*v1y)+(ii*v2y)
            t1z=az+(i*v1z)+(ii*v2z)
                        
                        
            t2x=bx+(i*v1x)-(pe*v1x)+(ii*v2x)
            t2y=by+(i*v1y)+(ii*v2y)
            t2z=bz+(i*v1z)+(ii*v2z)
                        
                        
            rt1x=(round(t1x, 6))
            rt1y=(round(t1y, 6))
            rt1z=(round(t1z, 6))
        
            rt2x=(round(t2x, 6))
            rt2y=(round(t2y, 6))
            rt2z=(round(t2z, 6))

            print string.rjust("C",0), string.rjust(`rt1x`, 12), string.rjust(`rt1y`, 12), string.rjust(`rt1z`, 12)
            print string.rjust("C",0), string.rjust(`rt2x`, 12), string.rjust(`rt2y`, 12), string.rjust(`rt2z`, 12)

##The odd numbered rows also require movement but by the same distance (both rows 5 & 6 require a leftward shift by one v1x
##distance) thus a second loop is created to handle the proper scaled leftward movements of these unit cells

        
        
    
for iii in range(m-1):
	tbotx= lax+((iii*v1x))
	rtbotx=(round(tbotx, 6))
	print string.rjust("C",0), string.rjust(`rtbotx`, 12), string.rjust(`lay`, 12), string.rjust(`laz`, 12)
        
##^This will create bottom zig-zag edge below row n=0 carbons to complete the benzene rings

if (n % 2 == 0):
	for i4 in range(m-1):
		ttopx= v1x+(i4*v1x)
    		ttopy= (n)*v2y
    		ttopz= lbz
    
    		rtopx= (round(ttopx, 6))
    		rtopy= (round(ttopy, 6))
    		rtopz= (round(ttopz, 6))
    
    		print string.rjust("C",0), string.rjust(`rtopx`, 12), string.rjust(`rtopy`, 12), string.rjust(`rtopz`, 12)
elif (n % 2 == 1):
	for i4 in range(m-1):
		ttopx= 1.3+(i4*v1x)
		ttopy= n*v2y
		ttopz= lbz

		rtopx= (round(ttopx, 6))
    		rtopy= (round(ttopy, 6))
    		rtopz= (round(ttopz, 6))
    
    		print string.rjust("C",0), string.rjust(`rtopx`, 12), string.rjust(`rtopy`, 12), string.rjust(`rtopz`, 12)
##^^This generates the carbons for the top zig-zag edge of the sheet

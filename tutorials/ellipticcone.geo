#
## An elliptic cone
#

# An elliptic cone, given by the point c = (c x , c y , c z ) at the base of the cone along the main axis,
# the vectors v and w of the long and short axis of the ellipse, respectively,
# the height of the cone, h, and ratio of base long axis length to top long axis length, r:
# ellipticcone (c x , c y , c z ; v x , v y , v z ; w x , w y , w z; h; r)

# Note: The elliptic cone has to be truncated by planes similar to a cone or an elliptic cylinder.
# When r =1, the truncated elliptic cone becomes an elliptic cylinder.
# When r tends to zero, the truncated elliptic cone tends to a full elliptic cone.
# However, when r = 0, the top part becomes a point(tip) and meshing fails!

algebraic3d

solid cutcone = ellipticcone ( 0, 0, 0; 5, 0, 0; 0, 2, 0; 5; 0.5)
	and plane (0, 0, 0; 0, 0, -1)
	and plane (0, 0, 5; 0, 0, 1);

tlo cutcone;

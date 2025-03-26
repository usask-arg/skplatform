.. _rotationmatrices:

Rotation Matrices
=================
Rotation matrices are used extensively in ``skplatform`` and although the subject is conceptually straight-forward it is very
easy to get confused with the subject.

Intrinsic and Extrinsic
-----------------------
The purpose of a rotation matrix is to rotate a vector to produce a new vector. This immediately opens up two common options:

    #. **Extrinsic**: hold the coordinate axes constant and rotate the given vector. The final vector is expressed in terms of the fixed axes.
    #. **Intrinsic**: hold the vector constant in the original coordinate system and rotate the coordinate axes to a new set of axes. The final vector is the original vector expressed in the new rotated corodinate system.

`Wikipedia Davenport Chained Rotations <https://en.wikipedia.org/wiki/Davenport_chained_rotations>`_ and `Wikipedia Euler Angles <https://en.wikipedia.org/wiki/Euler_angles>`_ describe the two rotations in much more detail.

Both types of rotation have application within ``skplatform``.

For example, specifying instrument look vectors in the :ref:`icf` are easily visualized using extrinsic rotations while orientian of the :ref:`pcf` using  yaw,
pitch and roll is best suited for intrinsic rotations.

Extrinsic Rotation Matrix
----------------------------------
Rotation matrices are used to (i) rotate vectors within one coordinate system or (ii) to transform the coordinates of a
vector in one coordinate system to coordinates in a rotated coordinate system.

Extrinsic rotations execute all the rotations in a sequence of rotations around the axes in the starting coordinate system. This starting
coordinate system is a *global* system that does not rotate.

The standard approach for a 3-axis cartesion coordinate system is to consider the effect of rotation around each of the
three axes upon each of the primary unit vectors: i.e. rotate the :math:`\hat{x}`, :math:`\hat{y}` and :math:`\hat{z}` unit vectors around the :math:`\hat{z}` axis
and use this to express the new, primed unit vectors in terms of the original unit vector.

For example, rotation around the :math:`\hat{z}` axis,

..  math::

    \hat{x}' &= &\cos\theta\:\hat{x}  &+ \sin{\theta}\:\hat{y} &+ 0\:\hat{z} \\
    \hat{y}' &= &-\sin\theta\:\hat{x} &+ \cos{\theta}\:\hat{y} &+ 0\:\hat{z} \\
    \hat{z}' &= &0\:\hat{x}           &+ 0\:\hat{y} &+ 1\:\hat{z} \\

We can then immediately write the extrinsic rotation matrix, :math:`R_z`, for rotation around the Z axis as,

..  math::

    R_{z} =  \begin{bmatrix}
                \cos\theta   & -\sin{\theta} & 0 \\
                \sin\theta &  \cos{\theta} & 0 \\
                 0          &  0            & 1
                \end{bmatrix}

and similary for Y and X gives,

..  math::

    R_{y} =  \begin{bmatrix}
                 \cos\theta & 0  & \sin{\theta} \\
                 0          & 1  & 0            \\
                -\sin\theta  & 0 & \cos{\theta}
                \end{bmatrix}

..  math::

    R_{x} =  \begin{bmatrix}
                1           & 0             & 0 \\
                0           & \cos\theta     & -\sin{\theta} \\
                0           & \sin\theta &  \cos{\theta}
                \end{bmatrix}


This rotation matrix will rotate any unit vector around the Z axis by angle :math:`\theta` and present the new unit vector
in terms of the original unit vectors.

For example: extrinsic rotation of the :math:`\hat{x}` unit vector is given by:

..  math::

        \hat{x}' = R_{z}\:\hat{x}

or more explicitly,

..  math::

    \hat{x}' =  \begin{bmatrix}
                \cos\theta   & -\sin{\theta} & 0 \\
                \sin\theta &  \cos{\theta} & 0 \\
                 0          &  0            & 1
                \end{bmatrix}
                \begin{bmatrix}
                1 \\
                0 \\
                0
                \end{bmatrix}

which upon evaluation gives

..  math::

    \hat{x}' =  \cos\theta\:\hat{x} + \sin{\theta}\:\hat{y} + 0\:\hat{z}

Intrinsic Rotation Matrix
-------------------------
Intrinsic rotations execute rotations in a sequence of rotations around the axes of a coordinate system after the axes have been rotated by
the previous rotation matrix. The rotating coordinate system is a *local* system that has a new orientation after application
of each rotation matrix. This rotating *local* system is easy to visualize and describe, especially for platforms like
aircraft, gondola, and spacecraft,  where yaw, pitch and roll are often used.

Theorem
^^^^^^^
Consider a sequence of rotations, *zyx*, applied to a vector :math:`\bar{v}`, where we first rotate around the Z axis,
then the Y axis and finally the X axis.

First consider the extrinsic case. In the extrinsic rotation the final vector, :math:`\bar{v}'_{ext}`, is given by the equation,

..  math::

    \bar{v}'_{ext} = \boldsymbol{R_x}\:\boldsymbol{R_y}\:\boldsymbol{R_z}\:\bar{v}

where we note that the order of rotations must be expressed mathematically from right to left so :math:`R_z` is first appled to :math:`\bar{v}`
then :math:`R_y` and finally :math:`R_x` and all rotations are around the fixed, starting, *global* coordinate system.

Now consider the intrinsic case. This is trickier to write down as the axes of rotation are changing between rotations, but
a proof is given below that proves the intrinsic rotation is equivalent to an extrinsic rotation with the order of rotation matrices reversed:

..  math::

    \boxed{\bar{v}'_{int} = \boldsymbol{R_z}\:\boldsymbol{R_y}\:\boldsymbol{R_x}\:\bar{v}}

Proof
^^^^^
`Wikipedia <http://en.wikipedia.org/wiki/Euler_angles#Conversion_between_intrinsic_and_extrinsic_rotations>`_ states

    Any extrinsic rotation is equivalent to an intrinsic rotation by the same angles but with inverted order of
    elemental rotations, and vice-versa. For instance, the intrinsic rotations x−y′−z′′ by angles α,β,γ are equivalent
    to the extrinsic rotations z−y−x by angles γ,β,α.

The rotation matrix, :math:`R`, can be considered to be a coordinate transform that expresses the unit vector :math:`\hat{x}'` in
terms of the unprimed unit vectors, :math:`\hat{x}`,  :math:`\hat{y}` and  :math:`\hat{z}`. Conversely, the inverse of the rotation matrix, :math:`R^{-1}`,
expresses the unprimed unit vectors, :math:`\hat{x}` in terms of the primed unit vectors :math:`\hat{x}'`.

Intrinsic rotations rely upon rotation around axes that are rotated along with the target vector, :math:`\bar{v}`. This presents an issue as
simply applying a rotation matrix to an unprimed vector:

..  math::

    \bar{v}'_g = \boldsymbol{R_z} \: \bar{v}_g

creates a new rotated vector, :math:`\bar{v}'_g`, where we haved used the *g* subscript to indicate the vector is expressed in terms of the original unprimed, *global*, unit vectors.
It is not appropriate to leave the rotated vector expressed in *global* unprimed components if we are going to perform a second (or third)
rotation around a *local*, primed axis. The *global* rotated vector must be written in terms of the new *local* unit vectors. This is achieved with a coordinate transform using the
inverse of the rotation matrix,

..  math::

    \bar{v}'_l = \boldsymbol{R_z^{-1}} \: \bar{v}'_g

Having transformed the primed vector to the local, rotated coordinate system we can now apply the second rotation, :math:`R_y` to produce
a double primed rotated vector,

..  math::

    \bar{v}''_l = \boldsymbol{R_y} \: \bar{v}'_l

This second,double primed, rotated vector is expressed in terms of the first set of *local* coordinates. To perform the third rotation
around the thrid axis we must transform the second, double primed  from the first *local* system to the second *local* system,

..  math::

    \bar{v}''_{ll} = \boldsymbol{R_y^{-1}} \: \bar{v}''_l

We have now expressed the second rotation vector in the second *local* coordinate system and we can now apply the third and final rotation,

..  math::

    \bar{v}'''_{ll} =\boldsymbol{R_x} \: \bar{v}''_{ll}

We now have the fully rotated vector but it is expressed in the second *local* system and we wish to get it back to the
original unprimed *global* system.  This is straight forward as we can easily transform the vector from the second, *local* system
to the first *local* system using the second rotation matrix,

..  math::

    \bar{v}'''_{l} =\boldsymbol{R_y} \: \bar{v}'''_{ll}

and from the first *local* system to the *global* system using the first rotation matrix,

..  math::

     \bar{v}'''_{g} =\boldsymbol{R_z} \: \bar{v}'''_{l}

Putting all the above equations together we get

..  math::

    \bar{v}'''_{g} = \boldsymbol{R_z} \: \boldsymbol{R_y} \: \boldsymbol{R_x} \: \boldsymbol{R_y^{-1}} \: \boldsymbol{R_y} \: \boldsymbol{R_z^{-1}} \: \boldsymbol{R_z}\: \bar{v}_{g}

or more succinctly

..  math::

    \boxed{ \bar{v}'''_{g} = \boldsymbol{R_z} \: \boldsymbol{R_y} \: \boldsymbol{R_x} \: \bar{v}_{g} }

And we have the proof of the standard, intrinsic rotation formula that rotation around the axes in the order z, y, x is implemented
by applying the (extrinsic like) rotation matrices in reverse order (from right to left), i.e. :math:`R_x` is applied first, :math:`R_y` second and
:math:`R_z` third.






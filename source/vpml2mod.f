!-----------------------------------------------------------------------
! This module contains the subroutines that allocate, setup and update all
! required information and arrays to attenuate the E and B fields using Vay's
! PML (perfectly matched layers) algorithm. It is based on his paper 
! JCP 165, 511 (2000).
! The subroutines are based off of the 2D subroutines from OSIRIS, located in
! os-emf-vpml.f90. A difference in this implementation is that the PML 
! boundaries to not "surround" the outside of the simulation box, but rather 
! overlap with it. Each PML boundary has its inner most cell x row or y collumn 
! copied to from the box, and all other cells are copied to the box (as opposed
! to just one row or collumn in OSIRIS). For a visual explaination of the 
! boundaries, see the map UPIC_vpml_maps.ods. That map is derived from the
! OSIRIS subroutines.
! The subroutines setup_vpml, update_e_vpml, and update_b_vpml, respectivly, are 
! called from the main new_bbeps2.f program and all other subroutines are
! called internally from this module.
! Originaly written by M.D. Meyers in June 2015.
! Last update July 25, 2015.
!
!
      module vpml2

      implicit none
      public :: setup_vpml
      public :: update_e_vpml, update_b_vpml
      private :: setup_att_coeff
      private :: update_interface_x, update_interface_y
      private :: update_e_x, update_e_y
      private :: update_b_x, update_b_y

      contains

         subroutine setup_vpml(vpmlx_e_up,vpmlx_e_low,vpmly_e_up, & 
                    vpmly_e_low,vpmlx_b_up,vpmlx_b_low,vpmly_b_up,vpmly_b_low,&
                    coef_e_low,coef_e_up,coef_b_low,coef_b_up,vpml_exyz, &
                    vpml_bxyz,vpml,xtrue,ytrue,nx,ny,nxe,nye,dt,ci)
! This subroutine allocates the vpml field arrays for each boundary
! It also calls the subroutine that calculates the attenuation coefficients
! in the pml layers.
! variables:
! vpmlx_e_up: E field in the right PML boundary.
! vpmlx_e_low: E field in the left pml boundary.
! vpmly_e_up: E field in the top PML boundary.
! vpmly_e_low: E field in the bottom pml boundary.
! coef_e_up(1:2,:,:): E field attenuation coefficients for the right (1)
! and top (2).
! coef_e_low(1:2,:,:): E field attenuation coefficients for the left (1)
! and bottom (2).
! vpml_exyz: This is a copy of the E field array inside the box.
! The same is true for all arrays with a "b" instead of an "e".
! vpml = the thickness of the pml boundary
! xtrue = (0,1) if x pml is (no,yes) present
! ytrue = (0,1) if y pml is (no,yes) present
! nx = physical x size of the simulation box partition
! ny = physical y size of the simulation box partition
! nxe = x size of the simulation box counting guard cells
! nye = y size of the simulation box couting guard cells
! dt = time step size.
! ci = inverse speed of light
         implicit none
         real, dimension(:,:,:), pointer :: vpmlx_e_up, vpmlx_e_low
         real, dimension(:,:,:), pointer :: vpmly_e_up, vpmly_e_low
         real, dimension(:,:,:), pointer :: vpmlx_b_up, vpmlx_b_low
         real, dimension(:,:,:), pointer :: vpmly_b_up, vpmly_b_low
         real, dimension(:,:,:), pointer :: coef_e_low, coef_e_up
         real, dimension(:,:,:), pointer :: coef_b_low, coef_b_up
         real, dimension(:,:,:), pointer :: vpml_exyz, vpml_bxyz
         integer :: vpml, xtrue, ytrue, nx, ny, nxe, nye
         real :: dt, ci

        ! unlike in OSIRIS, we make the PML overlap with the box

         if(xtrue==1) then ! if x pml

           allocate(vpmlx_e_low(4,1:vpml,1:ny))
           allocate(vpmlx_e_up(4,nx-vpml+1:nx,1:ny))

           allocate(vpmlx_b_low(4,1:vpml,1:ny))
           allocate(vpmlx_b_up(4,nx-vpml+1:nx,1:ny))

           vpmlx_e_up = 0.0; vpmlx_e_low = 0.0
           vpmlx_b_up = 0.0; vpmlx_b_low = 0.0

         endif

         if(xtrue==1.and.ytrue==1) then !if both x and y pml

           allocate(vpmly_e_low(4,vpml-1:nx-vpml+2,1:vpml))
           allocate(vpmly_e_up(4,vpml-1:nx-vpml+2,ny-vpml+1:ny))

           allocate(vpmly_b_low(4,vpml-1:nx-vpml+2,1:vpml))
           allocate(vpmly_b_up(4,vpml-1:nx-vpml+2,ny-vpml+1:ny))

           vpmly_e_up = 0.0; vpmly_e_low = 0.0
           vpmly_b_up = 0.0; vpmly_b_low = 0.0

         elseif(xtrue==0.and.ytrue==1) then !if only y pml

           allocate(vpmly_e_low(4,1:nx,1:vpml))
           allocate(vpmly_e_up(4,1:nx,ny-vpml+1:ny))

           allocate(vpmly_b_low(4,1:nx,1:vpml))
           allocate(vpmly_b_up(4,1:nx,ny-vpml+1:ny))

           vpmly_e_up = 0.0; vpmly_e_low = 0.0
           vpmly_b_up = 0.0; vpmly_b_low = 0.0

         endif

! allocate coefficient arrays
         allocate(coef_e_low(2,3,vpml),coef_e_up(2,3,vpml))
         allocate(coef_b_low(2,3,vpml),coef_b_up(2,3,vpml))

! allocate buffer box field arrays
         allocate(vpml_exyz(3,nxe,nye),vpml_bxyz(3,nxe,nye))

         vpml_exyz = 0.0; vpml_bxyz = 0.0

! calculate the attenuation coefficient arrays
         call setup_att_coeff(coef_e_low,coef_e_up,coef_b_low,coef_b_up,dt,ci)

         end subroutine setup_vpml

         subroutine setup_att_coeff(coef_e_low,coef_e_up,coef_b_low,coef_b_up,&
                    dt,ci)
! This subroutine calculates the pml field attenuation coefficients from 
! eq (27) of Vay's paper: JCP 165, 511 (2000)
! It is very similar to the OSIRIS routine with the same name.
! The first index in the coefficient arrays is (1,2) for (x,y) side.
! The suffex "low" is for the left or bottom bounday.
! The suffex "up" is for the right or top boundary.
! coef_e_low(1:2,:,:) = attenuation coefficients for E field inside 
! left (1) bottom (2) pml boundary.
! coef_e_up(1:2,:,:) = attenuation coefficients for E field inside 
! right (1) top (2) pml boundary.
! dt = time step size.
! ci = inverse speed of light
! thickness = pml thickness
! tj1, tj12, tj = amplitude attenuation factor between grid points 3/2 --> 1,
! 1 --> 1/2, 1/2 --> 0
! beta2e, beta2b = coefficients from Vay 2000, eq. (27)
! s1, s2 = attenuation exponents.
         implicit none
         real, dimension(:,:,:), pointer :: coef_e_low, coef_e_up
         real, dimension(:,:,:), pointer :: coef_b_low, coef_b_up
         real :: dt, ci
!local variables
         integer :: thickness, i, j
         real :: p_att_nexp = 4.6
         real :: tj, tj12, tj1, beta1e, beta1b
         real :: dtdx, dtdxdif  
         real :: beta2e, beta2b
         real :: s1, s2
         real :: dx

! The thickness of the PML (same for all boundaries).
         thickness = size(coef_e_low,3)

!Conversion factor between UPIC --> OSIRIS cell size
         dx = ci

! This gives optimal absorption according to eq (7.61) in Taflove
         s1 = - (0.8 * ( p_att_nexp + 1 )) / 2.0
         s2 = 1.0 / (thickness-3)
         dtdx = dt/dx

! The coefs are the same for the x and y boundaries. We could probably save
! a little memory in the future by exploiting this.
         do i=1, 2

	    ! (dx-dt)/(dx+dt) = (1-dt/dx)/(1+dt/dx)
	    dtdxdif = (1.0 - dtdx)/(1.0 + dtdx)
         
           do j=0, thickness-1

! Attenuation coefficients for LOWER boundaries [Eqs (38-40)] with \delta x = dx:
! tj = exp[-2 * (abs(j)/5)^n]
	  
! Check field positionings to obtain this (B is the outside point):
! j =  -0.5 => interface
! j =  -1.5 => one cell inside wall
! eq. 40 of vay 2000
! the loop starts from the gaurd cell at the inner boundary
	         tj   = exp( s1 * ( s2 * ( j - 0.5 ) )**p_att_nexp )
	         tj12 = exp( s1 * ( s2 * ( j ) )**p_att_nexp )
	         tj1  = exp( s1 * ( s2 * ( j + 0.5 ) )**p_att_nexp )              
 
! beta coefficients
! eq. 27 of vay 2000
	         beta1e = dtdx*( 1.0 + dtdxdif * (1.0 - tj12) )
	         beta2e = dtdx
	  
	         beta1b = dtdx*( 1.0 + dtdxdif * (1.0 - tj1 ) )
	         beta2b = dtdx

! coeffs. from E and B in eq. 27 of vay 2000
	         coef_e_low(i,1,thickness-j) = 1.0 - beta1e + tj12 * beta2e
	         coef_e_low(i,2,thickness-j) = tj * beta1e
	         coef_e_low(i,3,thickness-j) = beta2e
	
	         coef_b_low(i,1,thickness-j) = 1.0 - beta1b + tj1 * beta2b
	         coef_b_low(i,2,thickness-j) = tj12 * beta1b
	         coef_b_low(i,3,thickness-j) = beta2b

!Attenuation coefficients for UPPER boundaries
! E is the outside point
! j =  0 => interface
! j =  1 => one cell inside wall
	         tj   = exp( s1 * ( s2 *  j )**p_att_nexp )
	         tj12 = exp( s1 * ( s2 * (j + 0.5) )**p_att_nexp )
	         tj1  = exp( s1 * ( s2 * (j + 1.0) )**p_att_nexp )

! beta coefficients
	         beta1e = dtdx*( 1.0 + dtdxdif * (1.0 - tj12) )
	         beta1b = dtdx*( 1.0 + dtdxdif * (1.0 - tj1 ) )
 
	         coef_e_up(i,1,j+1) = 1.0 - beta1e + tj12 * dtdx
	         coef_e_up(i,2,j+1) = dtdx
	         coef_e_up(i,3,j+1) = tj * beta1e

	         coef_b_up(i,1,j+1) = 1.0 - beta1b + tj1 * dtdx
	         coef_b_up(i,2,j+1) = dtdx
	         coef_b_up(i,3,j+1) = tj12 * beta1b         

	       enddo
         enddo
  
         end subroutine setup_att_coeff

         subroutine update_e_vpml(vpmlx_e_up,vpmlx_e_low,vpmly_e_up, &
                    vpmly_e_low,vpmlx_b_up,vpmlx_b_low,vpmly_b_up,vpmly_b_low,&
                    b,coef_e_low,coef_e_up,xtrue,ytrue,dt,ci,inorder)
! this subroutine first communicates the B field between the pml layers
! and the box, and also between the x and y pml layers, if needed
! It then updates the E field in the PML layers.
! variables:
! vpmlx_e_low, vpmlx_e_up = E field inside (left,right) pml boundary.
! vpmly_e_low, vpmly_e_up = E field inside (bottom,top) pml boundary.
! Same as above for "b" instead of "e".
! b = box B field array
! coef_e_low(1:2,:,:) = attenuation coefficients for E field inside 
! left (1) bottom (2) pml boundary.
! coef_e_up(1:2,:,:) = attenuation coefficients for E field inside 
! right (1) top (2) pml boundary.
! xtrue = (0,1) if x pml is (no,yes) present
! ytrue = (0,1) if y pml is (no,yes) present
! inorder = interpolation order in the main code.
! dt = time step size.
! ci = inverse speed of light
! bxyz = pointer to 2-D box field array that's updated by the PML.
! side = (0,1) if (left or bottom,right or top) side.
! height = height of the box (not counting guard cells).
! box_width = the width of the box (not counting guard cells).
! thickness = pml thickness.
! dx = conversion factor between the OSIRIS and UPIC grid size.
         implicit none
         real, dimension(:,:,:), pointer :: vpmlx_e_low, vpmlx_e_up
         real, dimension(:,:,:), pointer :: vpmly_e_low, vpmly_e_up
         real, dimension(:,:,:), pointer :: vpmlx_b_low, vpmlx_b_up
         real, dimension(:,:,:), pointer :: vpmly_b_low, vpmly_b_up
         real, dimension(:,:,:), pointer :: coef_e_low, coef_e_up
         real, dimension(:,:,:), pointer :: coef_b_low, coef_b_up
         real, dimension(:,:,:), pointer :: b
         integer :: xtrue, ytrue, inorder
         real :: dt, ci
!local data
         real, dimension(:,:,:), pointer :: bxyz
         integer :: side !LHS or bottom = 0, RHS or top = 1
         integer :: height, box_width, thickness
         real :: dx

! Get the physical box height (y size) and box width (x size)
! Make the pointers to the 2D box field arrays point to the correct
! place in the actual box field buffer array
         if(inorder==1) then
           height = size(b,3) - 1
           box_width = size(b,2) - 2
           bxyz => b(:,1:box_width,1:height)
         elseif(inorder==2) then
           height = size(b,3) - 3
           box_width = size(b,2) - 4
           bxyz => b(:,2:box_width+1,2:height+1)
         elseif(inorder==3) then
           height = size(b,3) - 5
           box_width = size(b,2) - 6
           bxyz => b(:,3:box_width+2,3:height+2)
         endif
         
! get the thickness of the PML, and conversion factor between 
! UPIC --> OSIRIS cell size
         thickness = size(coef_e_low,3)
         dx = ci                  

! Update the boundaries between the parts of the x and y pml
! layers that border each other, and the pml borders with the box. 

         if(xtrue==1) then
      
           side = 0 ! left side
           call update_interface_x(vpmlx_b_low,vpmly_b_low,vpmly_b_up,bxyz, &
                box_width,height,thickness,ytrue,side)

           side = 1 ! right side
           call update_interface_x(vpmlx_b_up,vpmly_b_low,vpmly_b_up,bxyz, &
                box_width,height,thickness,ytrue,side)

         endif 

         if(ytrue==1) then

           side = 0 !bottom
           call update_interface_y(vpmly_b_low,bxyz,box_width,height, &
                thickness,xtrue,side)

           side = 1 !top
           call update_interface_y(vpmly_b_up,bxyz,box_width,height, &
                thickness,xtrue,side)

         endif

         ! Update the fields in the walls
         if(xtrue==1) then 

           side = 0 ! left side
           call update_e_x(vpmlx_e_low,vpmlx_b_low,coef_e_low,coef_e_up, &
                box_width,height,thickness,ytrue,side,dt,dx)

           side = 1 ! right side
           call update_e_x(vpmlx_e_up,vpmlx_b_up,coef_e_low,coef_e_up, &
                box_width,height,thickness,ytrue,side,dt,dx)

        endif

        if(ytrue==1) then

           side = 0 ! bottom
           call update_e_y(vpmly_e_low,vpmly_b_low,coef_e_low,box_width, &
                height,thickness,xtrue,side,dt,dx)

           side = 1 ! top
           call update_e_y(vpmly_e_up,vpmly_b_up,coef_e_up,box_width, &
                height,thickness,xtrue,side,dt,dx)

         endif

         end subroutine update_e_vpml

         subroutine update_interface_x(vpmlx_fld,vpmly_fld_low,vpmly_fld_up, &
                    box_fld,box_width,height,thickness,ytrue,side)
! This subroutine updates the interfaces between the box and x pml layers,
! and also between the x and y pml layers if needed.
! It is different from the OSIRIS subroutine in that it copies the fields from
! all pml cells, except the inner most collumn, to the box.
! It copies from the box to the inner most collumn of the pml only, liks OSIRIS.
! It also copies one cell from the x to y pml, and one cell from the y to x pml.
! variables:
! vpmlx_fld = field with being updated between interfaces.
! vpmly_fld_low, vpmly_fld_up = the lower or upper ypml field.
! box_fld = the field inside the box.
! box_width = the width of the box (not counting guard cells).
! height = height of the box (not counting guard cells).
! thickness = pml thickness.
! ytrue = (0,1) if y pml is (no,yes) present
! side = (0,1) if (left,right) side.
! xstart = x index of the left most cell being copied from the x pml to the box.
! xstop = x index of the right most cell being copied from the x pml to the box.
! box_to_pml = the x index of the column being copied from the box to the x pml.
! vpml_to_box = the x index of the inner most column being copied from
! the pml to the box.
! vpmlx_out = the x index of the cell whose field is being copied
! from the x pml (to some other array).
! vpmlx_in = the x index of the cell being copied to the x pml
! same as above for vpmly_out and vpmly_in
         implicit none
         real, dimension(:,:,:), pointer :: vpmlx_fld
         real, dimension(:,:,:), pointer :: vpmly_fld_low, vpmly_fld_up 
         real, dimension(:,:,:), pointer :: box_fld
         integer :: box_width, height, thickness, ytrue, side
!local variables
         integer :: xstart, xstop, box_to_vpml, vpml_to_box
         integer :: vpmlx_out, vpmlx_in, vpmly_out, vpmly_in
         integer :: box_edge
         integer :: ix, iy, iyg

! get the x loop range and cells to copy to/from the box for the correct x side
         if(side==0) then !left
           box_to_vpml = thickness
           vpml_to_box = thickness - 1
           box_edge = 1
           xstart = 1; xstop = vpml_to_box
         else !right
           box_to_vpml = box_width - thickness + 1
           vpml_to_box = box_to_vpml + 1
           box_edge = box_width
           xstart = vpml_to_box; xstop = box_width
         endif

! we always copy all the "interior" PML cells to the box. if ypmls are present, 
! then we need to copy to/from those too.
         vpmlx_out = vpml_to_box; vpmlx_in = box_to_vpml
         vpmly_out = box_to_vpml; vpmly_in = vpml_to_box

! now do the side of the box. start at the bottom. copy box field to vpml field
         do iy = 1, height
           vpmlx_fld(1:3,box_to_vpml,iy) = box_fld(1:3,box_to_vpml,iy)
           vpmlx_fld(4,box_to_vpml,iy) = 0.0
         enddo

! copy the fields to the box from the PML
         do ix = xstart, xstop
           do iy = 1, height
             box_fld(1:2,ix,iy) = vpmlx_fld(1:2,ix,iy)
             box_fld(3,ix,iy) = vpmlx_fld(3,ix,iy) + vpmlx_fld(4,ix,iy)
           enddo                    
         enddo

!         do iy = 1, height         
!           box_fld(:,box_edge,iy) = 0.0
!         enddo

! copy all field components to/from x pml from/to upper and lower y pml
         if(ytrue==1) then
         
           do iy = 1, thickness - 1 
             ! bottom
	         vpmlx_fld(:,vpmlx_in,iy) = vpmly_fld_low(:,vpmly_out,iy)
	         vpmly_fld_low(:,vpmly_in,iy) = vpmlx_fld(:,vpmlx_out,iy)
	         ! top, iyg is global y coordinate
             iyg = iy + height - thickness + 1 
	         vpmlx_fld(:,vpmlx_in,iyg) = vpmly_fld_up(:,vpmly_out,iyg)
	         vpmly_fld_up(:,vpmly_in,iyg) = vpmlx_fld(:,vpmlx_out,iyg)
           enddo

         endif

         end subroutine update_interface_x

         subroutine update_interface_y(vpmly_fld,box_fld,box_width,height, &
                    thickness,xtrue,side)
! This subroutine updates the interfaces between the box and y pml layers.
! It is different from the OSIRIS subroutine in that it copies the fields from
! all pml cells, except the inner most collumn, to the box.
! It copies from the box to the inner most collumn of the pml only, like OSIRIS.
! variables:
! vpmly_fld = field being updated between interfaces.
! box_fld = the field inside the box.
! box_width = the width of the box (not counting guard cells).
! height = height of the box (not counting guard cells).
! thickness = pml thickness.
! xtrue = (0,1) if x pml is (no,yes) present
! side = (0,1) if (bottom,top) side.
! xstart = x index of the left most pml column being updated.
! xstop = x index of the right most column being updated.
! ystart = y index of the lower most cell being copied from the y pml to the box.
! ystop = y index of the upper most cell being copied from the y pml to the box.
! box_to_pml = the y index of the row being copied from the box to the y pml.
! vpml_to_box = the y index of the inner most row being copied from
! the pml to the box.
         implicit none
         real, dimension(:,:,:), pointer :: vpmly_fld 
         real, dimension(:,:,:), pointer :: box_fld
         integer :: box_width, height, thickness, xtrue, side
!local variables
         integer :: xstart, xstop, ystart, ystop
         integer :: box_to_vpml, vpml_to_box
         integer :: ix, iy

! get the x loop range. same for both y pmls. loop over the whole y pml
         if(xtrue==0) then
           xstart = 1; xstop = box_width
         else
           xstart = thickness; xstop = box_width - thickness + 1
         endif

! get the y pml loop ranges for the correct y side. loop over whole y pml
         if(side==0) then !top
           box_to_vpml = thickness
           vpml_to_box = box_to_vpml - 1
           ystart = 1; ystop = vpml_to_box
         else !bottom
           box_to_vpml = height - (thickness - 1)
           vpml_to_box = box_to_vpml + 1
           ystart = vpml_to_box; ystop = height
         endif

! copy the pml fields to/from the box
! start at the left. copy box field to vpml field first. one column only
         do ix = xstart, xstop
           vpmly_fld(1:3,ix,box_to_vpml) = box_fld(1:3,ix,box_to_vpml)
           vpmly_fld(4,ix,box_to_vpml) = 0.0
         enddo

!copy vpml field to box field.
         do ix = xstart, xstop
           do iy = ystart, ystop
             box_fld(1:2,ix,iy) = vpmly_fld(1:2,ix,iy)
             box_fld(3,ix,iy) = vpmly_fld(3,ix,iy) + vpmly_fld(4,ix,iy)
           enddo      
         enddo 

         end subroutine update_interface_y

         subroutine update_e_x(vpmlx_e,vpmlx_b,coef_e_low,coef_e_up,box_width,&
                    height,thickness,ytrue,side,dt,dx)
! This subroutine updates the E field in the x pml layers according to
! UPIC_vpml_maps.ods.
! Note that this includes the corners (overlap between x and y pmls), if needed. 
! variables:
! vpmlx_e = E field in the pml being updated.
! vpmlx_b = B field for pml in which E is being updated.
! coef_e_low = attenuation coefficients for E in the left pml.
! coef_e_up = attenuation coefficients for E in the right pml.
! box_width = the width of the actual box (not counting guard cells).
! height = height of the box (not counting guard cells).
! thickness = pml thickness.
! ytrue = (0,1) if y pml is (no,yes) present
! side = (0,1) if (left,right) side.
! dt = time step size.
! dx = grid size (in OSIRIS units)
! coef_e = assigned to coef_e_low or coef_e_up depending on which side is
! being updated.
! coef_e_cnr_low = attenuation coefficients for the bottom corners.
! coef_e_cnr_up = attenuation coefficients for the upper corners.
! coef_e_cnr = assigned to coef_e_cnr_low or coef_e_cnr_up depending on which
! corner is being updated.
! dtdx = coeffecient for the finite difference solver.
! xstart = x index of the left most pml column being updated.
! xstart1 = x index of the left most invariant column being updated.
! xstop = x index of the right most column being updated.
! ystart_low, ystop_low = lower and upper y index of bottom corner regions.
! ystart_mid, ystop_mid = lower and upper y index of x pmls not in a corner.
! ystart_up, ystop_up = lower and upper y index of top corner regions.
! ystart = y index of the bottom of the invariant component loop range.
! ystop = y index of the top of the invariant component loop range.
         implicit none
         real, dimension(:,:,:), pointer :: vpmlx_e, vpmlx_b
         real, dimension(:,:,:), pointer :: coef_e_low, coef_e_up
         integer :: box_width, height, thickness, ytrue, side
         real :: dt, dx
!local vars
         real, dimension(:,:), pointer :: coef_e
         real, dimension(:,:), pointer :: coef_e_cnr_low, coef_e_cnr_up
         real, dimension(:), pointer :: coef_e_cnr
         real :: dtdx
         integer :: xstart, xstart1, xstop
         integer :: starty_low, stopy_low, starty_mid, stopy_mid
         integer :: starty_up, stopy_up, starty, stopy
         integer :: ix, iy, coef_i, coef_cnr_i

! coef. for finite difference solver
         dtdx = dt/dx

! get the coefs and pml range for the correct x side
         if(side==0) then !left
           coef_e => coef_e_low(1,:,:)
           xstart = 1; xstart1 = xstart + 1
           xstop = thickness - 1
         else !right
           coef_e => coef_e_up(1,:,:)
           xstart = box_width - thickness + 2; xstart1 = xstart
           xstop = box_width
         endif
         
! if there are y pmls, then we need to worry about the corners
! the loop ranges are always the same for both x sides
         if(ytrue==1) then
           starty = 1; stopy = height
           starty_low = 2; stopy_low = thickness - 1
           starty_mid = thickness; stopy_mid = height - thickness + 1
           starty_up = height - thickness + 2; stopy_up = height
         else
           starty = 1; stopy = height
           starty_low = 1; stopy_low = 0
           starty_up = 1; stopy_up = 0
           starty_mid = starty + 1; stopy_mid = stopy
         endif

! coefs for the y sides for the corners
         coef_e_cnr_low => coef_e_low(2,:,:)
         coef_e_cnr_up => coef_e_up(2,:,:)

! update the fields in the pml.
! loop through the regions for E in the UPIC_vpml_maps.ods
         do ix = xstart, xstop

! do the lower corner
           do iy = starty_low, stopy_low
             coef_cnr_i = iy
             coef_e_cnr => coef_e_cnr_low(:,coef_cnr_i)            
             !XY component => perp. solver (attenuation in corner)             
             vpmlx_e(1,ix,iy) = coef_e_cnr(1) * vpmlx_e(1,ix,iy) + &
                                coef_e_cnr(2) * ( vpmlx_b(3,ix,iy) +  &
                                vpmlx_b(4,ix,iy)) - coef_e_cnr(3) * &
                                ( vpmlx_b(3,ix,iy-1) + vpmlx_b(4,ix,iy-1) )
             !ZY component => perp. solver (attenuation in corner)
	         vpmlx_e(4,ix,iy) = coef_e_cnr(1) * vpmlx_e(4,ix,iy) - &
                                coef_e_cnr(2) * vpmlx_b(1,ix,iy) + &
                                coef_e_cnr(3) * vpmlx_b(1,ix,iy-1)
           enddo

! do the middle region that's not in a corner
           do iy = starty_mid , stopy_mid

             !XY component => standard Yee solver (no attenuton) 
             vpmlx_e(1,ix,iy) = vpmlx_e(1,ix,iy) + dtdx * ( vpmlx_b(3,ix,iy) &
                                + vpmlx_b(4,ix,iy)- vpmlx_b(3,ix,iy-1) &
                                - vpmlx_b(4,ix,iy-1) )
             !ZY component => standard Yee solver (no attenuton)
	         vpmlx_e(4,ix,iy) = vpmlx_e(4,ix,iy) - dtdx * ( vpmlx_b(1,ix,iy) &
	                            - vpmlx_b(1,ix,iy-1) )
           enddo

! do the upper corner
           do iy = starty_up, stopy_up
             coef_cnr_i = iy + thickness - height - 1
             coef_e_cnr => coef_e_cnr_up(:,coef_cnr_i)

             !XY component => perp. solver (attenuation in corner)             
             vpmlx_e(1,ix,iy) = coef_e_cnr(1) * vpmlx_e(1,ix,iy) + &
                                coef_e_cnr(2) * ( vpmlx_b(3,ix,iy) + &
                                vpmlx_b(4,ix,iy)) - coef_e_cnr(3) * &
                                ( vpmlx_b(3,ix,iy-1) + vpmlx_b(4,ix,iy-1) )
             !ZY component => perp. solver (attenuation in corner)                                
	         vpmlx_e(4,ix,iy) = coef_e_cnr(1) * vpmlx_e(4,ix,iy) - &
                                coef_e_cnr(2) * vpmlx_b(1,ix,iy) + &
                                coef_e_cnr(3) * vpmlx_b(1,ix,iy-1)
           enddo
         enddo

! Now do the invariant components, which are propagating only in x.
! Loop over the whole side.
         do ix = xstart1, xstop
           coef_i = ix + side*(thickness - box_width - 1)
           do iy = starty, stopy
           
             ! YX component => vpml wall solver (uniform attenuation in wall)
             vpmlx_e(2,ix,iy) = coef_e(1,coef_i) * vpmlx_e(2,ix,iy) - &
                                coef_e(2,coef_i) * ( vpmlx_b(3,ix,iy) + &
                                vpmlx_b(4,ix,iy) ) + coef_e(3,coef_i) * &
                                ( vpmlx_b(3,ix-1,iy) + vpmlx_b(4,ix-1,iy) )
             ! ZX component => vpml wall solver (uniform attenuation in wall)
             vpmlx_e(3,ix,iy) = coef_e(1,coef_i) * vpmlx_e(3,ix,iy) + &
                                coef_e(2,coef_i) * vpmlx_b(2,ix,iy) - &
                                coef_e(3,coef_i) * vpmlx_b(2,ix-1,iy)
           enddo
         enddo

         end subroutine update_e_x

         subroutine update_e_y(vpml_e,vpml_b,coef_e_in,box_width,height, &
                    thickness,xtrue,side,dt,dx)
! This subroutine updates the E field in the y pml layers according to
! UPIC_vpml_maps.ods.
! variables:
! vpmlx_e = E field in the pml being updated.
! vpmlx_b = B field for pml in which E is being updated.
! coef_e_in = passed in argument in coef_e_upper or coef_e_lower.
! height = height of a partition (not counting guard cells).
! thickness = pml thickness.
! xtrue = (0,1) if x pml is (no,yes) present
! side = (0,1) if (bottom,top) side.
! dt = time step size.
! dx = grid size (in OSIRIS units)
! coef_e = points to coef_e_in(2,:,:).
! dtdx = coeffecient for the finite difference solver.
! xstart = x index of the left most pml column being updated.
! xstart1 = (2,xstart) if (no,yes) x pml is present. This may be different from
! xstart because of the array indexes (they can't be less than 1).
! xstop = x index of the right most column being updated.
! ystart, ystop = lower and upper y index of the y pml region. Depends on xtrue.
! ystart1 = (2,ystart) if side = (0,1) . This may be different from ystart because
! of the array indexes (they can't be less than 1).
         implicit none
         real, dimension(:,:,:), pointer :: vpml_e, vpml_b
         real, dimension(:,:,:), pointer :: coef_e_in
         integer :: box_width, height, thickness, xtrue, side
         real :: dt, dx
! local variables
         real, dimension(:,:), pointer :: coef_e
         real :: dtdx
         integer :: xstart, xstart1, xstop, ystart, ystart1, ystop
         integer :: ix, iy, coef_i

! coef. for finite difference solver
         dtdx = dt/dx
         
! get the x loop range. same for both y pmls
         if(xtrue==0) then
           xstart = 1; xstart1 = xstart + 1
           xstop = box_width
         else
           xstart = thickness; xstart1 = xstart
           xstop = box_width - thickness + 1
         endif
         
! get the y pml loop ranges for the correct y side
         if(side==0) then
           ystart = 1; ystart1 = ystart + 1
           ystop = thickness - 1
         else
           ystart = height - thickness + 2; ystart1 = ystart
           ystop = height
         endif

! get the coefs for the y side
         coef_e => coef_e_in(2,:,:)

! loop through the regions for E in the UPIC_vpml_maps.ods
! The components propagating in y are attenuated
         do ix = xstart, xstop
           do iy = ystart1, ystop
             coef_i = iy + side*(thickness - height - 1)
             
             ! XY component => vpml wall solver (uniform attenuation in wall)
             vpml_e(1,ix,iy) = coef_e(1,coef_i) * vpml_e(1,ix,iy) + &
                               coef_e(2,coef_i) * ( vpml_b(3,ix,iy) + &
                               vpml_b(4,ix,iy) ) - coef_e(3,coef_i) * &
                               ( vpml_b(3,ix,iy-1) + vpml_b(4,ix,iy-1) )
             ! ZY component => vpml wall solver (uniform attenuation in wall)
             vpml_e(4,ix,iy) = coef_e(1,coef_i) * vpml_e(4,ix,iy) - &
                               coef_e(2,coef_i) * vpml_b(1,ix,iy) + &
                               coef_e(3,coef_i) * vpml_b(1,ix,iy-1)
           enddo
         enddo

! These need a different loop because of the index shift.
! the components propagating in x are not attenuated.
         do ix = xstart1, xstop
           do iy = ystart, ystop
           
             ! YX component => standard Yee solver (no attenuton)
             vpml_e(2,ix,iy) = vpml_e(2,ix,iy)- dtdx*( vpml_b(3,ix,iy) + &
                               vpml_b(4,ix,iy) - vpml_b(3,ix-1,iy) - &
                               vpml_b(4,ix-1,iy) )
             ! ZX component => standard Yee solver (no attenuton)
             vpml_e(3,ix,iy) = vpml_e(3,ix,iy) + dtdx*( vpml_b(2,ix,iy) - &
                               vpml_b(2,ix-1,iy) )
           enddo
         enddo

         end subroutine update_e_y

         subroutine update_b_vpml(vpmlx_b_up,vpmlx_b_low,vpmly_b_up, &
                    vpmly_b_low,vpmlx_e_up,vpmlx_e_low,vpmly_e_up,vpmly_e_low,&
                    f,coef_b_low,coef_b_up,xtrue,ytrue,dt,ci,inorder)
! this subroutine first communicates the E field between the pml layers
! and the box, and also between the x and y pml layers, if needed
! It then updates the B field in the PML layers.
! variables:
! vpmlx_b_low, vpmlx_b_up = B field inside (left,right) pml boundary.
! vpmly_b_low, vpmly_b_up = B field inside (bottom,top) pml boundary.
! Same as above for "e" instead of "b".
! f = box E field array
! coef_b_low(1:2,:,:) = attenuation coefficients for B field inside 
! left (1) bottom (2) pml boundary.
! coef_b_up(1:2,:,:) = attenuation coefficients for B field inside 
! right (1) top (2) pml boundary.
! xtrue = (0,1) if x pml is (no,yes) present
! ytrue = (0,1) if y pml is (no,yes) present
! inorder = interpolation order in the main code.
! dt = time step size.
! ci = inverse speed of light
! fxyz = pointer to 2-D box field array that's updated by the PML.
! side = (0,1) if (left or bottom,right or top) side.
! height = height of the box (not counting guard cells).
! box_width = the width of the box (not counting guard cells).
! thickness = pml thickness.
! dx = conversion factor between the OSIRIS and UPIC grid size.
         implicit none
         real, dimension(:,:,:), pointer :: vpmlx_b_low, vpmlx_b_up
         real, dimension(:,:,:), pointer :: vpmly_b_low, vpmly_b_up
         real, dimension(:,:,:), pointer :: vpmlx_e_low, vpmlx_e_up
         real, dimension(:,:,:), pointer :: vpmly_e_low, vpmly_e_up
         real, dimension(:,:,:), pointer :: coef_b_low, coef_b_up
         real, dimension(:,:,:), pointer :: coef_e_low, coef_e_up
         real, dimension(:,:,:), pointer :: f
         integer :: xtrue, ytrue, inorder
         real :: dt, ci
         !local variables
         real, dimension(:,:,:), pointer :: fxyz
         integer :: side !LHS or bottom = 0, RHS or top = 1
         integer :: height, box_width, thickness
         real :: dx

         if(inorder==1) then
           height = size(f,3) - 1
           box_width = size(f,2) - 2
           fxyz => f(:,1:box_width,1:height)
         elseif(inorder==2) then
           height = size(f,3) - 3
           box_width = size(f,2) - 4
           fxyz => f(:,2:box_width+1,2:height+1)
         elseif(inorder==3) then
           height = size(f,3) - 5
           box_width = size(f,2) - 6
           fxyz => f(:,3:box_width+2,3:height+2)
         endif

! get the thickness of the PML, and conversion factor between 
! UPIC --> OSIRIS cell size
         thickness = size(coef_b_low,3)
         dx = ci
         
! Convert the fields to OSIRIS units from UPIC units so we don't have to change
! the VPML equations.
         fxyz = fxyz*ci

! Update the boundaries between the parts of the x and y pml
! layers that border each other, and the pml borders with the box. 
         if(xtrue==1) then
      
           side = 0 ! left side
           call update_interface_x(vpmlx_e_low,vpmly_e_low,vpmly_e_up,fxyz, &
                box_width,height,thickness,ytrue,side)

           side = 1 ! right side
           call update_interface_x(vpmlx_e_up,vpmly_e_low,vpmly_e_up,fxyz, &
                box_width,height,thickness,ytrue,side)

         endif

         if(ytrue==1) then

           side = 0 !bottom
           call update_interface_y(vpmly_e_low,fxyz,box_width,height, &
                thickness,xtrue,side)

           side = 1 !top
           call update_interface_y(vpmly_e_up,fxyz,box_width,height, &
                thickness,xtrue,side)

         endif

! Update the fields in the walls
         if(xtrue==1) then

           side = 0 ! left side
           call update_b_x(vpmlx_b_low,vpmlx_e_low,coef_b_low,coef_b_up, &
                box_width,height,thickness,ytrue,side,dt,dx)

           side = 1 ! right side
           call update_b_x(vpmlx_b_up,vpmlx_e_up,coef_b_low,coef_b_up, &
                box_width,height,thickness,ytrue,side,dt,dx)

         endif

         if(ytrue==1) then

           side = 0 ! bottom
           call update_b_y(vpmly_b_low,vpmly_e_low,coef_b_low,box_width, &
                height,thickness,xtrue,side,dt,dx)

           side = 1 ! top
           call update_b_y(vpmly_b_up,vpmly_e_up,coef_b_up,box_width,height, &
                thickness,xtrue,side,dt,dx)

        endif

       ! Convert the fields back to UPIC units from OSIRIS units.
         fxyz = fxyz/ci

         end subroutine update_b_vpml

         subroutine update_b_x(vpmlx_b,vpmlx_e,coef_b_low,coef_b_up, &
                    box_width,height,thickness,ytrue,side,dt,dx)
! This subroutine updates the B field in the x pml layers according to
! UPIC_vpml_maps.ods.
! Note that this includes the corners (overlap between x and y pmls), if needed. 
! variables:
! vpmlx_b = B field in the pml being updated.
! vpmlx_e = E field for pml in which B is being updated.
! coef_b_low = attenuation coefficients for B in the left pml.
! coef_b_up = attenuation coefficients for B in the right pml.
! box_width = the width of the actual box (not counting guard cells).
! height = height of the box (not counting guard cells).
! thickness = pml thickness.
! ytrue = (0,1) if y pml is (no,yes) present
! side = (0,1) if (left,right) side.
! dt = time step size.
! dx = grid size (in OSIRIS units)
! coef_b = assigned to coef_b_low or coef_b_up depending on which side is
! being updated.
! coef_b_cnr_low = attenuation coefficients for the bottom corners.
! coef_b_cnr_up = attenuation coefficients for the upper corners.
! coef_b_cnr = assigned to coef_b_cnr_low or coef_e_cnr_up depending on which
! corner is being updated.
! dtdx = coeffecient for the finite difference solver.
! xstart = x index of the left most pml column being updated.
! xstop = x index of the right most column being updated.
! xstop1 = x index of the right most invariant column being updated.
! ystart_low, ystop_low = lower and upper y index of bottom corner regions.
! ystart_mid, ystop_mid = lower and upper y index of x pmls not in a corner.
! ystart_up, ystop_up = lower and upper y index of top corner regions.
! ystart = y index of the bottom of the invariant component loop range.
! ystop = y index of the top of the invariant component loop range.
         implicit none
         real, dimension(:,:,:), pointer :: vpmlx_b, vpmlx_e
         real, dimension(:,:,:), pointer :: coef_b_low, coef_b_up
         integer :: box_width, height, thickness, ytrue, side
         real :: dt, dx
         !local data
         real, dimension(:,:), pointer :: coef_b
         real, dimension(:,:), pointer :: coef_b_cnr_low, coef_b_cnr_up
         real, dimension(:), pointer :: coef_b_cnr
         real dtdx
         integer :: xstart, xstop, xstop1
         integer :: starty_low, stopy_low, starty_mid, stopy_mid
         integer :: starty_up, stopy_up, starty, stopy
         integer :: ix, iy, coef_i, coef_cnr_i

! coef. for finite difference solver
         dtdx = dt/dx

! get the coefs and pml range for the correct x side
         if(side==0) then
           coef_b => coef_b_low(1,:,:)
           xstart = 1
           xstop = thickness - 1; xstop1 = xstop
         else
           coef_b => coef_b_up(1,:,:)
           xstart = box_width - thickness + 2
           xstop = box_width; xstop1 = xstop - 1
           ! Note: fixed this is a bug. we really wanted xstop = box_width.
         endif

! if there are y pmls, then we need to worry about the corners
! the loop ranges are always the same for both x sides        
         if(ytrue==1) then
           starty = 1; stopy = height
           starty_low = 1; stopy_low = thickness - 1
           starty_mid = thickness; stopy_mid = height - thickness + 1
           starty_up = height - thickness + 2; stopy_up = height - 1
         else
           starty = 1; stopy = height
           starty_low = 1; stopy_low = 0
           starty_up = 1; stopy_up = 0
           starty_mid = starty; stopy_mid = stopy - 1
         endif

! coefs for the y sides for the corners
         coef_b_cnr_low => coef_b_low(2,:,:)
         coef_b_cnr_up => coef_b_up(2,:,:)

! update the fields in the pml.
! loop through the regions for B in the UPIC_vpml_maps.ods
         do ix = xstart, xstop

!do the lower corner
           do iy = starty_low, stopy_low 
             coef_cnr_i = iy + 1
             coef_b_cnr => coef_b_cnr_low(:,coef_cnr_i)
             !XY component => perp. solver (attenuation in corner) 
             vpmlx_b(1,ix,iy) = coef_b_cnr(1) * vpmlx_b(1,ix,iy) - &
                                 coef_b_cnr(2) * ( vpmlx_e(3,ix,iy+1) + &
                                 vpmlx_e(4,ix,iy+1) ) + coef_b_cnr(3) * &
                                 ( vpmlx_e(3,ix,iy) + vpmlx_e(4,ix,iy) ) 
             !ZY component => perp. solver (attenuation in corner)
	         vpmlx_b(4,ix,iy) = coef_b_cnr(1) * vpmlx_b(4,ix,iy) + &
	                            coef_b_cnr(2) * vpmlx_e(1,ix,iy+1) - &
	                            coef_b_cnr(3) * vpmlx_e(1,ix,iy)
           enddo

! do the middle region that's not in a corner
           do iy = starty_mid, stopy_mid 
             !XY component => standard Yee solver (no attenuton) 
             vpmlx_b(1,ix,iy) = vpmlx_b(1,ix,iy) - dtdx*( vpmlx_e(3,ix,iy+1) &
                                + vpmlx_e(4,ix,iy+1) - vpmlx_e(3,ix,iy) - &
                                vpmlx_e(4,ix,iy) )
             !ZY component => standard Yee solver (no attenuton)  
	         vpmlx_b(4,ix,iy) = vpmlx_b(4,ix,iy) + dtdx*(vpmlx_e(1,ix,iy+1) - &
                                vpmlx_e(1,ix,iy) )
           enddo

! do the upper corner
           do iy = starty_up, stopy_up
             coef_cnr_i = iy + thickness - (height + 1)
             coef_b_cnr => coef_b_cnr_up(:,coef_cnr_i)

             !XY component => perp. solver (attenuation in corner) 
             vpmlx_b(1,ix,iy) = coef_b_cnr(1) * vpmlx_b(1,ix,iy) - &
                                coef_b_cnr(2) * ( vpmlx_e(3,ix,iy+1) + &
                                vpmlx_e(4,ix,iy+1) ) + coef_b_cnr(3) * &
                                ( vpmlx_e(3,ix,iy) + vpmlx_e(4,ix,iy) ) 
             !ZY component => perp. solver (attenuation in corner) 
	         vpmlx_b(4,ix,iy) = coef_b_cnr(1) * vpmlx_b(4,ix,iy) + &
	                            coef_b_cnr(2) * vpmlx_e(1,ix,iy+1) - &
	                            coef_b_cnr(3) * vpmlx_e(1,ix,iy)
           enddo
         enddo

! Now do the invariant components, which are propagating only in x.
! Loop over the whole side.
         do ix = xstart, xstop1
           coef_i = ix + 1 + side*(thickness - box_width - 2)
           do iy = starty, stopy
     
             ! YX component => vpml wall solver (uniform attenuation in wall)
             vpmlx_b(2,ix,iy) = coef_b(1,coef_i) * vpmlx_b(2,ix,iy) + &
                                coef_b(2,coef_i) * ( vpmlx_e(3,ix+1,iy) + &
                                vpmlx_e(4,ix+1,iy) ) - coef_b(3,coef_i) * &
                                ( vpmlx_e(3,ix,iy) + vpmlx_e(4,ix,iy) )
             ! ZX component => vpml wall solver (uniform attenuation in wall)
             vpmlx_b(3,ix,iy) = coef_b(1,coef_i) * vpmlx_b(3,ix,iy) - &
                                coef_b(2,coef_i) * vpmlx_e(2,ix+1,iy) + &
                                coef_b(3,coef_i) * vpmlx_e(2,ix,iy) 
           enddo
         enddo

         end subroutine update_b_x

         subroutine update_b_y(vpml_b,vpml_e,coef_b_in,box_width,height, &
                    thickness,xtrue,side,dt,dx)
! This subroutine updates the B field in the y pml layers according to
! UPIC_vpml_maps.ods.
! variables:
! vpmlx_b = B field in the pml being updated.
! vpmlx_e = E field for pml in which B is being updated.
! coef_b_in = passed in argument in coef_b_upper or coef_b_lower.
! height = height of a partition (not counting guard cells).
! thickness = pml thickness.
! xtrue = (0,1) if x pml is (no,yes) present
! side = (0,1) if (bottom,top) side.
! dt = time step size.
! dx = grid size (in OSIRIS units)
! coef_b = points to coef_b_in(2,:,:).
! dtdx = coeffecient for the finite difference solver.
! xstart = x index of the left most pml column being updated.
! xstop = x index of the right most column being updated.
! xstop1 = (xstop,xstop+1) if (no,yes) x pml is present. This may be different from
! xstart because of the array indexes (they can't be less than 1).
! ystart, ystop = lower and upper y index of the y pml region. Depends on xtrue.
! ystart1 = (2,ystart) if side = (0,1) . This may be different from ystart because
! of the array indexes (they can't be less than 1).
! ystop1 = (ystop,ystop-1) if side = (0,1) . This may be different from ystop because
! of the array indexes (they can't be more than height).
         implicit none
         real, dimension(:,:,:), pointer :: vpml_b, vpml_e
         real, dimension(:,:,:), pointer :: coef_b_in
         integer :: box_width, height, thickness, xtrue, side
         real :: dt, dx
!local data
         real, dimension(:,:), pointer :: coef_b
         real :: dtdx
         integer :: xstart, xstop, xstop1, ystart, ystop, ystop1
         integer :: ix, iy, coef_i

! coef. for finite difference solver
         dtdx = dt/dx
         
! get the x loop range. same for both y pmls
         if(xtrue==0) then
           xstart = 1
           xstop = box_width; xstop1 = xstop - 1
         else
           xstart = thickness
           xstop = box_width - thickness + 1; xstop1 = xstop
         endif

! get the y pml loop ranges for the correct y side
         if(side==0) then !bottom
           ystart = 1
           ystop = thickness - 1; ystop1 = ystop
         else ! top
           ystart = height - thickness + 2
           ystop = height - 1; ystop1 = ystop - 1
         endif

! get the coefs for the y side
         coef_b => coef_b_in(2,:,:)

! loop through the regions for B in the UPIC_vpml_maps.ods
! The components propagating in y are attenuated
         do ix = xstart, xstop
           do iy = ystart, ystop1
             coef_i = iy + 1 + side*(thickness - height - 2)

             ! XY component => vpml wall solver (uniform attenuation in wall)
             vpml_b(1,ix,iy) = coef_b(1,coef_i) * vpml_b(1,ix,iy) - &
                               coef_b(2,coef_i) * ( vpml_e(3,ix,iy+1) + &
                               vpml_e(4,ix,iy+1) ) + coef_b(3,coef_i) * &
                               ( vpml_e(3,ix,iy) + vpml_e(4,ix,iy) )
             ! ZY component => vpml wall solver (uniform attenuation in wall)
             vpml_b(4,ix,iy) = coef_b(1,coef_i) * vpml_b(4,ix,iy) + &
                               coef_b(2,coef_i) * vpml_e(1,ix,iy+1) - &
                               coef_b(3,coef_i) * vpml_e(1,ix,iy)
           enddo
         enddo

! These need a different loop because of the index shift.
! the components propagating in x are not attenuated.
         do ix = xstart, xstop1
           do iy = ystart, ystop

             ! YX component => standard Yee solver (no attenuton)
             vpml_b(2,ix,iy) = vpml_b(2,ix,iy) + dtdx*( vpml_e(3,ix+1,iy) + &
                               vpml_e(4,ix+1,iy) - vpml_e(3,ix,iy) - &
                               vpml_e(4,ix,iy) )
             ! ZX component => standard Yee solver (no attenuton)
             vpml_b(3,ix,iy) = vpml_b(3,ix,iy) - dtdx*( vpml_e(2,ix+1,iy) - &
                               vpml_e(2,ix,iy) )
           enddo
         enddo

         end subroutine update_b_y

      end module vpml2 



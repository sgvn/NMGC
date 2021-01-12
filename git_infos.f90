!******************************************************************************
! MODULE: git_infos
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Automatically generated module (with Makefile.py) that contain
!! information on code version and branch name. 
!
!> @warning Do not modify this file manually !
!
!******************************************************************************
module git_infos
implicit none

character(len=40), parameter :: commit = '2fd5048eba7e368158ea7cefcb06fe27f2649dc1' !< commit ID when binary was compiled 
character(7), parameter :: branch = 'develop' !< branch name when compilation occured
character(len=80), parameter :: modifs = '/!\ There is non committed modifications'

end module git_infos
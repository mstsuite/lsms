! Created by F.Moitzi on 17.12.2022.

program test_mixing
   use PublicTypeDefinitionsModule, only: MixListStruct
   use MixingModule
   implicit none

   integer :: size = 3
   integer :: nt = 1
   integer :: nq(1) = 1

   real(kind=8), target, allocatable :: vector_old(:)
   real(kind=8), target, allocatable :: vector_new(:)
   real(kind=8), target, allocatable :: mesh(:)

   type(MixListStruct) :: mix_list

   allocate(vector_old(size), source = 0.0_8)
   allocate(vector_new(size), source = 0.0_8)
   allocate(mesh(size), source = 0.0_8)

   mix_list % vector_old => vector_old
   mix_list % vector_new => vector_new
   mix_list % mesh => mesh
   mix_list % size = size
   mix_list % rms = 0.0_8

   call initMixing(nt, nq, mix_list)

   call setBroydenMixing(1, 1, 0.02_8)

   vector_old = [0.1_8, 0.15_8, 0.2_8]
   vector_new = vector_old

   !-- 1.

   vector_new(1) = vector_new(1) + 0.01_8
   vector_new(2) = vector_new(2) + 0.02_8
   vector_new(3) = vector_new(3) - 0.01_8

   mix_list % rms = sqrt(sum(vector_new ** 2 - vector_old ** 2))

   print *, mix_list % rms

   call mixValues(mix_list)

   print '(3f20.8)', vector_old
   print '(3f20.8)', vector_new

   vector_old = vector_new

   mix_list % vector_old => vector_old
   mix_list % vector_new => vector_new
   mix_list % mesh => mesh
   mix_list % size = size

   !-- 2.

   vector_new(1) = vector_new(1) + 0.001_8
   vector_new(2) = vector_new(2) + 0.002_8
   vector_new(3) = vector_new(3) - 0.001_8

   mix_list % rms = sqrt(sum(vector_new ** 2 - vector_old ** 2))

   print *, mix_list % rms

   call mixValues(mix_list)

   print '(3f20.8)', vector_old
   print '(3f20.8)', vector_new

   vector_old = vector_new
   mix_list % vector_old => vector_old
   mix_list % vector_new => vector_new
   mix_list % mesh => mesh
   mix_list % size = size

   !-- 3.

   vector_new(1) = vector_new(1) + 0.0001_8
   vector_new(2) = vector_new(2) + 0.0002_8
   vector_new(3) = vector_new(3) - 0.0001_8

   mix_list % rms = sqrt(sum(vector_new ** 2 - vector_old ** 2))

   call mixValues(mix_list)

   print '(3f20.8)', vector_old
   print '(3f20.8)', vector_new

   vector_old = vector_new


   print *, "THE END"

   call endMixing()


end program test_mixing
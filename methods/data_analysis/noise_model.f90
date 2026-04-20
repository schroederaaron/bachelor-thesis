module noise_model
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use, intrinsic :: iso_fortran_env, only: int32, real64
  implicit none
  
  integer, parameter :: dp = real64
  integer, parameter :: ip = int32
  
  ! Derived type for sorted data
  type :: sorted_data_t
    real(dp), allocatable :: means_sorted(:)
    integer, allocatable :: original_indices(:)
    integer, allocatable :: n_residuals(:)
    real(dp), allocatable :: residuals_packed(:,:)  ! max_resid x n_genes
    integer :: max_resid_per_gene
    integer :: n_genes
  end type sorted_data_t
  
  ! Cache for family pools
  type :: family_cache_t
    real(dp), allocatable :: family_pools(:,:)  ! max_pool_size x n_families
    real(dp), allocatable :: orth_pools(:,:)
    integer, allocatable :: family_pool_sizes(:)
    integer, allocatable :: orth_pool_sizes(:)
    logical, allocatable :: is_cached(:)
  end type family_cache_t
  
contains

  ! ============================================================
  ! Helper functions for quicksort
  ! ============================================================
  pure function real_less(a, b) result(less)
    real(dp), intent(in) :: a, b
    logical :: less
    less = a < b
  end function real_less
  
  pure function real_greater(a, b) result(greater)
    real(dp), intent(in) :: a, b
    logical :: greater
    greater = a > b
  end function real_greater
  
  pure subroutine swap_int(a, b)
    integer(ip), intent(inout) :: a, b
    integer(ip) :: temp
    temp = a
    a = b
    b = temp
  end subroutine swap_int

  ! ============================================================
  ! Internal quicksort implementation for real arrays.
  ! Sorts indirectly using the permutation vector perm. Manual stack replaces recursion.
  ! ============================================================
  pure subroutine quicksort_real(array, perm, n, stack_left, stack_right)
    ! Real input array to sort
    real(dp), intent(in) :: array(:)
    ! Permutation vector that will be sorted
    integer(ip), intent(inout) :: perm(:)
    ! Size of the array
    integer(ip), intent(in) :: n
    ! Manual stack of left indices for quicksort recursion
    integer(ip), intent(inout) :: stack_left(:)
    ! Manual stack of right indices for quicksort recursion
    integer(ip), intent(inout) :: stack_right(:)

    integer(ip) :: left, right, i, j, top, pivot_idx
    real(dp) :: pivot_val

    top = 1
    stack_left(top) = 1
    stack_right(top) = n

    ! Iterative quicksort using explicit stack
    do while (top > 0)
      left = stack_left(top)
      right = stack_right(top)
      top = top - 1

      if (left >= right) cycle

      ! Select pivot and initialize pointers
      pivot_idx = (left + right) / 2
      pivot_val = array(perm(pivot_idx))
      i = left
      j = right

      ! Partitioning loop
      do
        do while (real_less(array(perm(i)), pivot_val))
          i = i + 1
        end do
        do while (real_greater(array(perm(j)), pivot_val))
          j = j - 1
        end do
        if (i <= j) then
          call swap_int(perm(i), perm(j))
          i = i + 1
          j = j - 1
        end if
        if (i > j) exit
      end do

      ! Push new ranges onto stack
      if (left < j) then
        top = top + 1
        stack_left(top) = left
        stack_right(top) = j
      end if
      if (i < right) then
        top = top + 1
        stack_left(top) = i
        stack_right(top) = right
      end if
    end do
  end subroutine quicksort_real

  ! ============================================================
  ! Random number utilities
  ! ============================================================
  subroutine init_random_seed()
    integer :: n
    integer, allocatable :: seed(:)
    
    call random_seed(size=n)
    allocate(seed(n))
    seed = 42
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine init_random_seed
  
  ! ============================================================
  ! Sort means and pack residuals for efficient access
  ! ============================================================
  subroutine prepare_sorted_data(means, replicates, n_samples, n_genes, &
                                 sorted_data, max_resid_per_gene)
    real(dp), intent(in) :: means(:)
    real(dp), intent(in) :: replicates(:,:)  ! n_samples x n_genes
    integer, intent(in) :: n_samples, n_genes
    type(sorted_data_t), intent(out) :: sorted_data
    integer, intent(out) :: max_resid_per_gene
    
    integer :: i, j, idx
    integer, allocatable :: order(:)
    integer, allocatable :: stack_left(:), stack_right(:)
    integer :: stack_size
    real(dp) :: gene_mean
    
    ! Store dimensions
    sorted_data%n_genes = n_genes
    
    ! Allocate order array
    allocate(order(n_genes))
    
    ! Initialize order with original indices
    do i = 1, n_genes
      order(i) = i
    end do
    
    ! Allocate stack for quicksort (maximum depth needed: 2 * log2(n) + some buffer)
    stack_size = 2 * int(log(real(n_genes, dp)) / log(2.0_dp)) + 10
    allocate(stack_left(stack_size), stack_right(stack_size))
    
    ! Sort order by means using quicksort
    call quicksort_real(means, order, n_genes, stack_left, stack_right)
    
    deallocate(stack_left, stack_right)
    
    ! Count residuals per gene and find maximum
    max_resid_per_gene = n_samples
    
    ! Allocate sorted arrays
    allocate(sorted_data%means_sorted(n_genes))
    allocate(sorted_data%original_indices(n_genes))
    allocate(sorted_data%n_residuals(n_genes))
    allocate(sorted_data%residuals_packed(max_resid_per_gene, n_genes))
    
    ! Fill sorted data
    do i = 1, n_genes
      idx = order(i)
      sorted_data%original_indices(i) = idx
      sorted_data%means_sorted(i) = means(idx)
      sorted_data%n_residuals(i) = n_samples
      
      !! Fill with "centered" residuals
      gene_mean = sum(replicates(:, idx)) / real(n_samples, dp)
      do j = 1, n_samples
        sorted_data%residuals_packed(j, i) = replicates(j, idx) - gene_mean
      end do
    end do
    
    sorted_data%max_resid_per_gene = max_resid_per_gene
    
    deallocate(order)
    
  end subroutine prepare_sorted_data

  ! ============================================================
  ! Binary search to find closest position in sorted array
  ! ============================================================
  function find_closest(target, means_sorted) result(pos)
    real(dp), intent(in) :: target
    real(dp), intent(in) :: means_sorted(:)
    integer :: pos, left, right, mid, n
    
    n = size(means_sorted)
    
    ! Handle edge cases
    if (n == 0) then
      pos = 0
      return
    end if
    
    if (target <= means_sorted(1)) then
      pos = 1
      return
    end if
    
    if (target >= means_sorted(n)) then
      pos = n
      return
    end if
    
    ! Binary search
    left = 1
    right = n
    
    do while (left <= right)
      mid = (left + right) / 2
      if (means_sorted(mid) < target) then
        left = mid + 1
      else if (means_sorted(mid) > target) then
        right = mid - 1
      else
        pos = mid
        return
      end if
    end do
    
    ! Find closest position
    if (left == 1) then
      pos = 1
    else if (left > n) then
      pos = n
    else
      if (abs(means_sorted(left) - target) < abs(means_sorted(left-1) - target)) then
        pos = left
      else
        pos = left - 1
      end if
    end if
    
  end function find_closest

  ! ============================================================
  ! Gather residuals adaptively - OPTIMIZED
  ! ============================================================
  subroutine gather_residuals_optimized(target_mean, sorted_data, &
                                      k_start, k_step, k_max, tau, &
                                      pooled_residuals, n_pooled, max_pool_size, &
                                      temp_pool)
    
    real(dp), intent(in) :: target_mean
    type(sorted_data_t), intent(in) :: sorted_data
    integer, intent(in) :: k_start, k_step, k_max
    real(dp), intent(in) :: tau
    real(dp), intent(out) :: pooled_residuals(:)
    integer, intent(out) :: n_pooled
    integer, intent(in) :: max_pool_size
    real(dp), intent(inout) :: temp_pool(:)  ! Pre-allocated workspace
    
    integer :: pos, left_cand, right_cand
    integer :: idx
    integer :: current_size, added_this_round, genes_added, offset
    real(dp) :: S_old, S_new, rel_change
    integer :: pool_size, n_resid
    
    ! Initialize
    n_pooled = 0
    current_size = 0
    
    ! Find closest position
    pos = find_closest(target_mean, sorted_data%means_sorted)
    if (pos == 0) return
    
    ! Initialize with closest gene
    current_size = min(sorted_data%n_residuals(pos), max_pool_size)
    pooled_residuals(1:current_size) = sorted_data%residuals_packed(1:current_size, pos)
    
    ! Initialize neighbors
    left_cand = pos - 1
    right_cand = pos + 1
    
    ! Initial expansion to reach k_start
    do while (current_size < k_start .and. (left_cand >= 1 .or. right_cand <= sorted_data%n_genes))
      ! Choose next nearest neighbor
      if (left_cand >= 1 .and. right_cand <= sorted_data%n_genes) then
        if (abs(sorted_data%means_sorted(left_cand) - target_mean) <= &
            abs(sorted_data%means_sorted(right_cand) - target_mean)) then
          idx = left_cand
          left_cand = left_cand - 1
        else
          idx = right_cand
          right_cand = right_cand + 1
        end if
      else if (left_cand >= 1) then
        idx = left_cand
        left_cand = left_cand - 1
      else
        idx = right_cand
        right_cand = right_cand + 1
      end if
      
      ! Add residuals from this gene
      pool_size = min(current_size + sorted_data%n_residuals(idx), max_pool_size)
      call add_residuals_to_pool(pooled_residuals, current_size, &
                                 sorted_data%residuals_packed(:, idx), &
                                 sorted_data%n_residuals(idx), pool_size)
      current_size = pool_size
    end do
    
    ! Check if we have enough residuals
    if (current_size < 10) return
    
    ! Compute initial mean absolute residual
    S_old = sum(abs(pooled_residuals(1:current_size))) / real(current_size, dp)
    
    if (S_old == 0.0_dp) then
      n_pooled = current_size
      return
    end if
    
    ! Adaptive growth
    do while (current_size < min(k_max, max_pool_size) .and. &
              (left_cand >= 1 .or. right_cand <= sorted_data%n_genes))
      
      added_this_round = 0
      genes_added = 0
      offset = 0  ! Track cumulative addition within this round
      
      ! Use temp_pool as workspace
      temp_pool(1:current_size) = pooled_residuals(1:current_size)
      
      ! Add neighbors until we have at least k_step new residuals
      do while (added_this_round < k_step .and. &
                (left_cand >= 1 .or. right_cand <= sorted_data%n_genes) .and. &
                current_size + offset < min(k_max, max_pool_size))
        
        ! Choose next nearest neighbor
        if (left_cand >= 1 .and. right_cand <= sorted_data%n_genes) then
          if (abs(sorted_data%means_sorted(left_cand) - target_mean) <= &
              abs(sorted_data%means_sorted(right_cand) - target_mean)) then
            idx = left_cand
            left_cand = left_cand - 1
          else
            idx = right_cand
            right_cand = right_cand + 1
          end if
        else if (left_cand >= 1) then
          idx = left_cand
          left_cand = left_cand - 1
        else
          idx = right_cand
          right_cand = right_cand + 1
        end if
        
        ! Get number of residuals for this gene
        n_resid = sorted_data%n_residuals(idx)
        
        ! Add residuals from this gene - use current_size + offset as starting position
        pool_size = min(current_size + offset + n_resid, max_pool_size)
        call add_residuals_to_pool(temp_pool, current_size + offset, &
                                   sorted_data%residuals_packed(:, idx), &
                                   n_resid, pool_size)
        
        ! Update counters
        offset = offset + n_resid
        added_this_round = added_this_round + n_resid
        genes_added = genes_added + 1
      end do
      
      if (genes_added == 0) exit
      
      ! Compute new mean absolute residual (total size = current_size + offset)
      S_new = sum(abs(temp_pool(1:current_size + offset))) / &
              real(current_size + offset, dp)
      rel_change = (S_new - S_old) / S_old
      
      ! Check stopping condition
      if (rel_change > tau) exit
      
      ! Accept expansion
      current_size = current_size + offset
      pooled_residuals(1:current_size) = temp_pool(1:current_size)
      S_old = S_new
    end do
    
    n_pooled = current_size
    
  end subroutine gather_residuals_optimized
  
  ! Helper subroutine to add residuals to pool
  subroutine add_residuals_to_pool(pool, current_size, residuals, n_resid, new_size)
    real(dp), intent(inout) :: pool(:)
    integer, intent(in) :: current_size
    real(dp), intent(in) :: residuals(:)
    integer, intent(in) :: n_resid, new_size
    
    integer :: n_to_copy
    
    n_to_copy = min(n_resid, new_size - current_size)
    if (n_to_copy > 0) then
      pool(current_size+1:current_size+n_to_copy) = residuals(1:n_to_copy)
    end if
  end subroutine add_residuals_to_pool

  ! ============================================================
  ! Monte Carlo p-value computation with fixed B
  ! ============================================================
  subroutine compute_pvalue_fixed(mu_c, r_c, resid_c, n_c, &
                                          mu_h, r_h, resid_h, n_h, &
                                          obs, B, norm_method, p)
    
    use, intrinsic :: iso_fortran_env, only: int64, real64
    implicit none
    
    real(real64), intent(in) :: mu_c, mu_h, obs
    integer, intent(in) :: r_c, r_h, B, n_c, n_h, norm_method
    real(real64), intent(in) :: resid_c(:), resid_h(:)
    real(real64), intent(out) :: p
    
    ! Local variables
    real(real64), allocatable :: null_dists(:)
    real(real64), allocatable :: random_c(:,:), random_h(:,:)
    integer :: count_ge, i, j
    real(real64) :: log2_factor, obs_abs
    real(real64) :: eta_c, eta_h, x_c, x_h, null_dist
    logical :: use_log_transform
    integer :: idx_c, idx_h
    
    ! ============================================================
    ! Precompute constants
    ! ============================================================
    use_log_transform = (norm_method /= 0)
    if (use_log_transform) then
        ! Precompute 1/log(2) instead of dividing each time
        log2_factor = 1.0_real64 / log(2.0_real64)
    end if
    
    ! Precompute absolute observed value (avoid recomputing in loop)
    obs_abs = abs(obs)
    
    ! ============================================================
    ! Allocate arrays
    ! ============================================================
    allocate(null_dists(B))
    allocate(random_c(r_c, B), random_h(r_h, B))
    
    ! ============================================================
    ! OPTIMIZATION 1: Batch random number generation
    ! Instead of B * (r_c + r_h) calls, we use just 2 calls
    ! ============================================================
    call random_number(random_c)  ! Generates r_c * B random numbers
    call random_number(random_h)  ! Generates r_h * B random numbers
    
    count_ge = 0
    
    ! ============================================================
    ! OPTIMIZATION 2: Main loop with manual optimizations
    ! - Direct index calculation without temporary arrays
    ! - Fused operations to reduce temporary variables
    ! - Branchless clamping where possible
    ! ============================================================
    do i = 1, B
        
        ! Compute eta_c (mean of sampled cancer residuals)
        eta_c = 0.0_real64
        do j = 1, r_c
            ! Convert random number to index in [1, n_c]
            ! Using multiplication and int() is faster than floor()
            idx_c = int(random_c(j, i) * real(n_c, real64)) + 1
            
            ! Branchless clamping using min/max intrinsics
            ! These are optimized by gfortran to single instructions
            idx_c = min(max(idx_c, 1), n_c)
            
            eta_c = eta_c + resid_c(idx_c)
        end do
        eta_c = eta_c / real(r_c, real64)
        
        ! Compute eta_h (mean of sampled healthy residuals)
        eta_h = 0.0_real64
        do j = 1, r_h
            idx_h = int(random_h(j, i) * real(n_h, real64)) + 1
            idx_h = min(max(idx_h, 1), n_h)
            eta_h = eta_h + resid_h(idx_h)
        end do
        eta_h = eta_h / real(r_h, real64)
        
        ! ============================================================
        ! OPTIMIZATION 3: Fused arithmetic operations
        ! Combine operations to reduce temporary variables
        ! ============================================================
        x_c = mu_c + eta_c
        x_h = mu_h + eta_h
        
        if (use_log_transform) then
            ! Use max() intrinsic - very fast in gfortran
            x_c = max(x_c, 0.0_real64) + 1.0_real64
            x_h = max(x_h, 0.0_real64) + 1.0_real64
            
            ! Use precomputed log2_factor instead of dividing by log(2)
            null_dist = abs((log(x_c) - log(x_h)) * log2_factor)
        else
            null_dist = abs(x_c - x_h)
        end if
        
        null_dists(i) = null_dist
        
        ! Count using integer comparison (faster than real comparison)
        if (null_dist >= obs_abs) count_ge = count_ge + 1
    end do
    
    ! ============================================================
    ! Compute final p-value using integer arithmetic
    ! ============================================================
    p = real(count_ge + 1, real64) / real(B + 1, real64)
    
    ! Cleanup
    deallocate(null_dists, random_c, random_h)
    
  end subroutine compute_pvalue_fixed

end module noise_model

  ! ============================================================
  ! Main wrapper for R interface - FULLY OPTIMIZED
  ! ============================================================
  subroutine compute_noise_pvalues_fortran( &
      cancer_means, cancer_replicates, cancer_n_genes, cancer_n_samples, &
      healthy_means, healthy_replicates, healthy_n_genes, healthy_n_samples, &
      obs_own, obs_fam, obs_orth, &
      family_means, ortholog_means, &
      valid_genes_own, valid_genes_fam, valid_genes_orth, &
      family_sizes, is_ortholog_sum, gene_to_family, &
      n_genes, n_families, norm_method, B, k_start, k_step, k_max, tau, &
      pvalues_own, pvalues_fam, pvalues_orth, n_success, &
      max_pool_size, neighborhood_size_own, neighborhood_size_fam, neighborhood_size_orth, neighborhood_size_cancer)
    
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use noise_model
    implicit none
    
    ! Input arrays
    real(c_double), intent(in) :: cancer_means(cancer_n_genes)
    real(c_double), intent(in) :: cancer_replicates(cancer_n_samples, cancer_n_genes)
    integer(c_int), intent(in) :: cancer_n_genes, cancer_n_samples
    
    real(c_double), intent(in) :: healthy_means(healthy_n_genes)
    real(c_double), intent(in) :: healthy_replicates(healthy_n_samples, healthy_n_genes)
    integer(c_int), intent(in) :: healthy_n_genes, healthy_n_samples
    
    real(c_double), intent(in) :: obs_own(n_genes), obs_fam(n_genes), obs_orth(n_genes)
    real(c_double), intent(in) :: family_means(n_families), ortholog_means(n_families)
    integer(c_int), intent(in) :: family_sizes(n_families)
    integer(c_int), intent(in) :: is_ortholog_sum
    integer(c_int), intent(in) :: gene_to_family(n_genes)
    integer(c_int), intent(in) :: n_genes, n_families
    integer(c_int), intent(in) :: norm_method, B, k_start, k_step, k_max
    real(c_double), intent(in) :: tau
    integer(c_int), intent(in) :: max_pool_size
    integer(c_int), intent(in) :: valid_genes_own(n_genes)
    integer(c_int), intent(in) :: valid_genes_fam(n_genes)
    integer(c_int), intent(in) :: valid_genes_orth(n_genes)
    
    ! Output
    real(c_double), intent(out) :: pvalues_own(n_genes), pvalues_fam(n_genes), pvalues_orth(n_genes)
    integer(c_int), intent(out) :: neighborhood_size_own(n_genes), neighborhood_size_fam(n_genes), neighborhood_size_orth(n_genes), neighborhood_size_cancer(n_genes)
    integer(c_int), intent(out) :: n_success
    
    ! Local variables
    integer :: i, fam_id, r_c, r_h, family_size, n_orth
    real(dp) :: mu_c, mu_h, family_mean, ortholog_mean
    real(dp) :: obs_own_val, obs_fam_val, obs_orth_val
    real(dp), allocatable :: cancer_resid_pool(:), healthy_resid_pool_own(:)
    real(dp), allocatable :: temp_pool(:)  ! Workspace
    integer :: cancer_resid_count, healthy_resid_count_own
    
    ! Sorted data structures
    type(sorted_data_t) :: cancer_sorted, healthy_sorted
    integer :: max_resid_per_gene
    
    ! Cache for family pools - PRE-COMPUTED
    type(family_cache_t) :: cache
    integer :: fam
    
    ! Initialize random seed
    call init_random_seed()
    
    ! Prepare sorted data
    call prepare_sorted_data(cancer_means, cancer_replicates, &
                            cancer_n_samples, cancer_n_genes, &
                            cancer_sorted, max_resid_per_gene)
    
    call prepare_sorted_data(healthy_means, healthy_replicates, &
                            healthy_n_samples, healthy_n_genes, &
                            healthy_sorted, max_resid_per_gene)
    
    ! Pre-compute family and ortholog pools
    allocate(cache%family_pools(max_pool_size, n_families))
    allocate(cache%orth_pools(max_pool_size, n_families))
    allocate(cache%family_pool_sizes(n_families))
    allocate(cache%orth_pool_sizes(n_families))
    allocate(cache%is_cached(n_families))
    cache%is_cached = .false.
    
    ! Pre-allocate workspace
    allocate(temp_pool(max_pool_size * 2))
    
    ! Pre-compute all family pools
    do fam = 1, n_families
      if (family_sizes(fam) > 0) then
        call gather_residuals_optimized(family_means(fam), healthy_sorted, &
                                        k_start, k_step, k_max, tau, &
                                        cache%family_pools(:, fam), &
                                        cache%family_pool_sizes(fam), &
                                        max_pool_size, temp_pool)
        
        call gather_residuals_optimized(ortholog_means(fam), healthy_sorted, &
                                        k_start, k_step, k_max, tau, &
                                        cache%orth_pools(:, fam), &
                                        cache%orth_pool_sizes(fam), &
                                        max_pool_size, temp_pool)
        cache%is_cached(fam) = .true.
      end if
    end do
    
    ! Process each gene
    n_success = 0
    r_c = cancer_n_samples
    r_h = healthy_n_samples
    n_orth = is_ortholog_sum
    pvalues_fam = -1.0_dp
    pvalues_orth = -1.0_dp
    pvalues_own = -1.0_dp
    
    neighborhood_size_cancer = -1
    neighborhood_size_fam = -1
    neighborhood_size_orth = -1
    neighborhood_size_own = -1
    
    ! Pre-allocate per-gene arrays (reused)
    allocate(cancer_resid_pool(max_pool_size))
    allocate(healthy_resid_pool_own(max_pool_size))
    
    do i = 1, n_genes

      ! Get values
      mu_c = cancer_means(i)
      mu_h = healthy_means(i)
      fam_id = gene_to_family(i)
      
      if (fam_id < 1 .or. fam_id > n_families) cycle
      if (.not. cache%is_cached(fam_id)) cycle
      
      family_mean = family_means(fam_id)
      ortholog_mean = ortholog_means(fam_id)
      family_size = family_sizes(fam_id)
      

      ! Gather pools (reusing pre-allocated arrays)
      call gather_residuals_optimized(mu_c, cancer_sorted, &
                                      k_start, k_step, k_max, tau, &
                                      cancer_resid_pool, cancer_resid_count, &
                                      max_pool_size, temp_pool)
      if (cancer_resid_count < 10) cycle
      
      call gather_residuals_optimized(mu_h, healthy_sorted, &
                                      k_start, k_step, k_max, tau, &
                                      healthy_resid_pool_own, healthy_resid_count_own, &
                                      max_pool_size, temp_pool)
      if (healthy_resid_count_own < 10) cycle

      
      ! Get observed differences
      obs_own_val = obs_own(i)
      obs_fam_val = obs_fam(i)
      obs_orth_val = obs_orth(i)
      
      if (obs_own_val /= obs_own_val .or. &
          obs_fam_val /= obs_fam_val .or. &
          obs_orth_val /= obs_orth_val) cycle
      
      ! Compute p-values with fixed B      
      if(valid_genes_own(i) == 1) then
        call compute_pvalue_fixed(mu_c, r_c, cancer_resid_pool(1:cancer_resid_count), &
                                cancer_resid_count, &
                                mu_h, r_h, healthy_resid_pool_own(1:healthy_resid_count_own), &
                                healthy_resid_count_own, &
                                obs_own_val, B, norm_method, &
                                pvalues_own(i))
        neighborhood_size_own(i) = healthy_resid_count_own
      end if
      
      if(valid_genes_fam(i) == 1) then
        call compute_pvalue_fixed(mu_c, r_c, cancer_resid_pool(1:cancer_resid_count), &
                                cancer_resid_count, &
                                family_mean, r_h, &
                                cache%family_pools(1:cache%family_pool_sizes(fam_id), fam_id), &
                                cache%family_pool_sizes(fam_id), &
                                obs_fam_val, B, norm_method, &
                                pvalues_fam(i))
        neighborhood_size_fam(i) = cache%family_pool_sizes(fam_id)
      end if
      
      if(valid_genes_orth(i) == 1) then
        call compute_pvalue_fixed(mu_c, r_c, cancer_resid_pool(1:cancer_resid_count), &
                                cancer_resid_count, &
                                ortholog_mean, r_h, &
                                cache%orth_pools(1:cache%orth_pool_sizes(fam_id), fam_id), &
                                cache%orth_pool_sizes(fam_id), &
                                obs_orth_val, B, norm_method, &
                                pvalues_orth(i))  
        neighborhood_size_orth(i) = cache%orth_pool_sizes(fam_id)
      end if

      neighborhood_size_cancer(i) = cancer_resid_count

      n_success = n_success + 1
    end do
    
    ! Clean up
    deallocate(cancer_resid_pool, healthy_resid_pool_own, temp_pool)
    deallocate(cache%family_pools, cache%orth_pools)
    deallocate(cache%family_pool_sizes, cache%orth_pool_sizes, cache%is_cached)
    
    if (allocated(cancer_sorted%means_sorted)) deallocate(cancer_sorted%means_sorted)
    if (allocated(cancer_sorted%original_indices)) deallocate(cancer_sorted%original_indices)
    if (allocated(cancer_sorted%n_residuals)) deallocate(cancer_sorted%n_residuals)
    if (allocated(cancer_sorted%residuals_packed)) deallocate(cancer_sorted%residuals_packed)
    
    if (allocated(healthy_sorted%means_sorted)) deallocate(healthy_sorted%means_sorted)
    if (allocated(healthy_sorted%original_indices)) deallocate(healthy_sorted%original_indices)
    if (allocated(healthy_sorted%n_residuals)) deallocate(healthy_sorted%n_residuals)
    if (allocated(healthy_sorted%residuals_packed)) deallocate(healthy_sorted%residuals_packed)
    
  end subroutine compute_noise_pvalues_fortran
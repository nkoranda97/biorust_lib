/// Parallel map: apply `$f` to each element of `$slice`, collecting into a Vec.
macro_rules! par_map {
    ($slice:expr, $f:expr) => {{
        #[cfg(feature = "parallel")]
        {
            use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
            $slice.par_iter().map($f).collect()
        }
        #[cfg(not(feature = "parallel"))]
        {
            $slice.iter().map($f).collect()
        }
    }};
}

/// Parallel fallible map: apply `$f` returning Result to each element, collecting into Result<Vec>.
macro_rules! par_try_map {
    ($slice:expr, $f:expr) => {{
        #[cfg(feature = "parallel")]
        {
            use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
            $slice.par_iter().map($f).collect()
        }
        #[cfg(not(feature = "parallel"))]
        {
            $slice.iter().map($f).collect()
        }
    }};
}

/// Parallel mutable for-each: apply `$f` to each element of `$slice` in place.
macro_rules! par_for_each_mut {
    ($slice:expr, $f:expr) => {{
        #[cfg(feature = "parallel")]
        {
            use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
            $slice.par_iter_mut().for_each($f);
        }
        #[cfg(not(feature = "parallel"))]
        {
            $slice.iter_mut().for_each($f);
        }
    }};
}

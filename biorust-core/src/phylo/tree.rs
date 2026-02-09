use crate::error::{BioError, BioResult};

use super::distance::DistanceMatrix;

#[derive(Debug, Clone)]
pub struct PhyloNode {
    pub label: Option<Box<str>>,
    pub branch_length: Option<f64>,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct PhyloTree {
    nodes: Vec<PhyloNode>,
    root: Option<usize>,
}

impl PhyloTree {
    pub fn root(&self) -> Option<usize> {
        self.root
    }

    pub fn node(&self, idx: usize) -> &PhyloNode {
        &self.nodes[idx]
    }

    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn num_leaves(&self) -> usize {
        self.nodes.iter().filter(|n| n.children.is_empty()).count()
    }

    pub fn leaves(&self) -> Vec<usize> {
        self.nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.children.is_empty())
            .map(|(i, _)| i)
            .collect()
    }

    pub fn leaf_labels(&self) -> Vec<String> {
        self.nodes
            .iter()
            .filter(|n| n.children.is_empty())
            .map(|n| n.label.as_deref().unwrap_or("").to_string())
            .collect()
    }

    pub fn nodes(&self) -> &[PhyloNode] {
        &self.nodes
    }
}

pub fn neighbor_joining(dist: &DistanceMatrix) -> BioResult<PhyloTree> {
    let n = dist.n();
    if n < 2 {
        return Err(BioError::TooFewSequences { n });
    }

    // Arena: n leaves + up to (n-2) internal nodes
    let max_nodes = 2 * n - 2;
    let mut nodes: Vec<PhyloNode> = Vec::with_capacity(max_nodes);

    // Create leaf nodes
    for label in dist.labels() {
        nodes.push(PhyloNode {
            label: Some(label.clone()),
            branch_length: None,
            parent: None,
            children: Vec::new(),
        });
    }

    // Working distance matrix: size (2n-1) x (2n-1) to accommodate new nodes
    let cap = 2 * n - 1;
    let mut d = vec![0.0f64; cap * cap];
    for i in 0..n {
        for j in 0..n {
            d[i * cap + j] = dist.get(i, j);
        }
    }

    // Active node indices
    let mut active: Vec<usize> = (0..n).collect();
    let mut next_node = n;

    while active.len() > 2 {
        let r = active.len();

        // Compute row sums
        let mut row_sum = vec![0.0f64; cap];
        for &i in &active {
            for &j in &active {
                row_sum[i] += d[i * cap + j];
            }
        }

        // Find minimum Q
        let mut min_q = f64::INFINITY;
        let mut min_i = 0;
        let mut min_j = 0;
        for (ai, &i) in active.iter().enumerate() {
            for &j in &active[(ai + 1)..] {
                let q = (r as f64 - 2.0) * d[i * cap + j] - row_sum[i] - row_sum[j];
                if q < min_q {
                    min_q = q;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        // Branch lengths to new node
        let dij = d[min_i * cap + min_j];
        let r_f = r as f64;
        let li = dij / 2.0 + (row_sum[min_i] - row_sum[min_j]) / (2.0 * (r_f - 2.0));
        let lj = dij - li;

        // Create new internal node
        let u = next_node;
        next_node += 1;
        nodes.push(PhyloNode {
            label: None,
            branch_length: None,
            parent: None,
            children: vec![min_i, min_j],
        });
        nodes[min_i].parent = Some(u);
        nodes[min_i].branch_length = Some(li);
        nodes[min_j].parent = Some(u);
        nodes[min_j].branch_length = Some(lj);

        // Update distances
        for &k in &active {
            if k == min_i || k == min_j {
                continue;
            }
            let duk = (d[min_i * cap + k] + d[min_j * cap + k] - dij) / 2.0;
            d[u * cap + k] = duk;
            d[k * cap + u] = duk;
        }
        d[u * cap + u] = 0.0;

        // Update active set: remove min_i, min_j; add u
        active.retain(|&x| x != min_i && x != min_j);
        active.push(u);
    }

    // Final two nodes: connect them
    debug_assert_eq!(active.len(), 2);
    let a = active[0];
    let b = active[1];
    let dab = d[a * cap + b];

    // For NJ (unrooted), we create one more internal node connecting the last two
    let u = next_node;
    nodes.push(PhyloNode {
        label: None,
        branch_length: None,
        parent: None,
        children: vec![a, b],
    });
    nodes[a].parent = Some(u);
    nodes[a].branch_length = Some(dab / 2.0);
    nodes[b].parent = Some(u);
    nodes[b].branch_length = Some(dab / 2.0);

    Ok(PhyloTree {
        nodes,
        root: None, // NJ is unrooted
    })
}

pub fn upgma(dist: &DistanceMatrix) -> BioResult<PhyloTree> {
    let n = dist.n();
    if n < 2 {
        return Err(BioError::TooFewSequences { n });
    }

    // Arena: n leaves + (n-1) internal nodes
    let max_nodes = 2 * n - 1;
    let mut nodes: Vec<PhyloNode> = Vec::with_capacity(max_nodes);

    for label in dist.labels() {
        nodes.push(PhyloNode {
            label: Some(label.clone()),
            branch_length: None,
            parent: None,
            children: Vec::new(),
        });
    }

    let cap = 2 * n - 1;
    let mut d = vec![0.0f64; cap * cap];
    for i in 0..n {
        for j in 0..n {
            d[i * cap + j] = dist.get(i, j);
        }
    }

    let mut active: Vec<usize> = (0..n).collect();
    let mut cluster_size = vec![1usize; cap];
    let mut heights = vec![0.0f64; cap];
    let mut next_node = n;

    while active.len() > 1 {
        // Find minimum distance pair
        let mut min_d = f64::INFINITY;
        let mut min_i = 0;
        let mut min_j = 0;
        for (ai, &i) in active.iter().enumerate() {
            for &j in &active[(ai + 1)..] {
                if d[i * cap + j] < min_d {
                    min_d = d[i * cap + j];
                    min_i = i;
                    min_j = j;
                }
            }
        }

        let u = next_node;
        next_node += 1;
        let h = min_d / 2.0;
        heights[u] = h;

        nodes.push(PhyloNode {
            label: None,
            branch_length: None,
            parent: None,
            children: vec![min_i, min_j],
        });
        nodes[min_i].parent = Some(u);
        nodes[min_i].branch_length = Some(h - heights[min_i]);
        nodes[min_j].parent = Some(u);
        nodes[min_j].branch_length = Some(h - heights[min_j]);

        let si = cluster_size[min_i];
        let sj = cluster_size[min_j];
        cluster_size[u] = si + sj;

        // Weighted average distances
        for &k in &active {
            if k == min_i || k == min_j {
                continue;
            }
            let duk = (d[min_i * cap + k] * si as f64 + d[min_j * cap + k] * sj as f64)
                / (si + sj) as f64;
            d[u * cap + k] = duk;
            d[k * cap + u] = duk;
        }
        d[u * cap + u] = 0.0;

        active.retain(|&x| x != min_i && x != min_j);
        active.push(u);
    }

    let root = active[0];
    Ok(PhyloTree {
        nodes,
        root: Some(root),
    })
}

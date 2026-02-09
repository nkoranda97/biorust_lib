use super::tree::PhyloTree;

pub fn to_newick(tree: &PhyloTree) -> String {
    // Find the serialization root: either the tree root (UPGMA) or the last node (NJ pseudo-root)
    let start = tree
        .root()
        .unwrap_or_else(|| tree.num_nodes().saturating_sub(1));

    let mut s = String::new();
    write_subtree(tree, start, &mut s);
    s.push(';');
    s
}

fn needs_quoting(label: &str) -> bool {
    label.chars().any(|ch| {
        ch.is_whitespace() || matches!(ch, ':' | ',' | '(' | ')' | ';' | '[' | ']' | '\'')
    })
}

fn write_label(out: &mut String, label: &str) {
    if label.is_empty() {
        return;
    }
    if needs_quoting(label) {
        out.push('\'');
        for ch in label.chars() {
            if ch == '\'' {
                out.push_str("''");
            } else {
                out.push(ch);
            }
        }
        out.push('\'');
    } else {
        out.push_str(label);
    }
}

fn write_subtree(tree: &PhyloTree, idx: usize, out: &mut String) {
    let node = tree.node(idx);

    if node.children.is_empty() {
        // Leaf node
        if let Some(ref label) = node.label {
            write_label(out, label);
        }
    } else {
        out.push('(');
        for (i, &child) in node.children.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            write_subtree(tree, child, out);
            if let Some(bl) = tree.node(child).branch_length {
                out.push(':');
                out.push_str(&format!("{:.6}", bl));
            }
        }
        out.push(')');
        if let Some(ref label) = node.label {
            write_label(out, label);
        }
    }
}

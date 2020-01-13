
ad_dir = "/Users/jeremys/work/opentargets/AD_finemap"
df = read_tsv(file.path(ad_dir, "network/ad_network.averages.all.tsv"))
candidates = read_tsv(file.path(ad_dir, "candidate_genes.tsv", col_names = "symbol"))

df$candidate = NA
df$candidate[df$gene %in% candidates$symbol] = df$gene[df$gene %in% candidates$symbol]
df$is_candidate = !is.na(df$candidate)
df = df %>% arrange(desc(rankingIte1000.node.median)) %>% mutate(index = row_number())
p1 = ggplot(mapping = aes(x=index, y=rankingIte1000.node.median, col=is_candidate, size=is_candidate)) +
  geom_point(data = df %>% filter(!is_candidate), alpha = 0.1) +
  geom_point(data = df %>% filter(is_candidate), alpha = 0.7) +
  scale_size_manual(values = c("FALSE" = 1, "TRUE" = 3)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "cornflowerblue")) +
  geom_text(data = df %>% filter(is_candidate), aes(label = candidate), nudge_x = 1500, size=2, col="black")
p1

p1.zoom = ggplot(mapping = aes(x=index, y=rankingIte1000.node.median, col=is_candidate, size=is_candidate)) +
  geom_point(data = df %>% filter(!is_candidate), alpha = 0.1) +
  geom_point(data = df %>% filter(is_candidate), alpha = 0.7) +
  scale_size_manual(values = c("FALSE" = 0.5, "TRUE" = 3)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "cornflowerblue")) +
  geom_text(data = df %>% filter(is_candidate), aes(label = candidate), nudge_x = 150, size=2, col="black") +
  coord_cartesian(xlim=c(0,2000), ylim=c(90, 100))
p1.zoom

pdf(file.path(ad_dir, "plots/network_candidate_genes.pdf"), width=7, height=6)
p1
p1.zoom
dev.off()

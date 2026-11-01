"""Load the JSON config and turn relative paths into absolute paths."""

import json
import os


PATH_KEYS = {
    "output_dir",
    "fasta_dir",
    "busco_root",
    "transcriptome_list",
    "masked_cds_fasta",
    "spliced_cds_fasta",
    "alignment_dir",
    "gff",
    "transcript_map",
    "species_map",
    "thresholds_file",
}

STAGES = [
    "transcriptome_qc",
    "busco_summary",
    "longest_isoform",
    "species_map",
    "reference_cds_check",
    "gff_metrics",
    "exon_qc",
    "msa_metrics",
    "compare_gff_msa",
    "filter",
    "report",
]

THRESHOLD_KEYS = [
    "n_exons_gt",
    "n_introns_lt",
    "max_exon_bp_gt",
    "max_intron_bp_lt",
    "avg_exon_bp_gt",
    "avg_intron_bp_lt",
]


class ConfigError(Exception):
    pass


def load_config(path):
    path = os.path.abspath(path)
    with open(path) as handle:
        cfg = json.load(handle)
    if not isinstance(cfg, dict):
        raise ConfigError("Config root must be a JSON object.")
    _resolve(cfg, os.path.dirname(path))
    cfg["_config_path"] = path
    if "output_dir" not in cfg:
        raise ConfigError("Config is missing output_dir.")
    if "alignment_suffix" not in cfg:
        raise ConfigError("Config is missing alignment_suffix.")
    if cfg.get("reference_row") != "last":
        raise ConfigError(
            'Config reference_row must be "last" '
            "(MarkerMiner writes the reference sequence as the last FASTA record)."
        )
    return cfg


def _resolve(node, base):
    if isinstance(node, dict):
        for key, value in list(node.items()):
            if key in PATH_KEYS and isinstance(value, str) and value:
                node[key] = _abs(base, value)
            else:
                _resolve(value, base)
    elif isinstance(node, list):
        for item in node:
            _resolve(item, base)


def _abs(base, path):
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(base, path))


def section(cfg, name):
    value = cfg.get(name) or {}
    if not isinstance(value, dict):
        raise ConfigError("Config section {0} must be an object.".format(name))
    return value


def enabled(cfg, name):
    return bool(section(cfg, name).get("enabled"))


def require(section_dict, key, label):
    if key not in section_dict or section_dict[key] in ("", None):
        raise ConfigError("{0} is missing {1}.".format(label, key))
    return section_dict[key]


def _check_threshold_keys(thresholds, label):
    if not isinstance(thresholds, dict):
        raise ConfigError("{0} must be a JSON object.".format(label))
    missing = [key for key in THRESHOLD_KEYS if key not in thresholds]
    if missing:
        raise ConfigError("{0} is missing: {1}".format(label, ", ".join(missing)))
    return thresholds


def resolve_thresholds(filter_section, required):
    """Cutoffs come from thresholds_file when it is set, otherwise from thresholds.

    The file is the JSON downloaded from report.html. A filter run needs one
    of the two. The report can be built before either exists.
    """
    path = filter_section.get("thresholds_file") or ""
    if path:
        if not os.path.isfile(path):
            raise ConfigError("filter.thresholds_file not found: {0}".format(path))
        with open(path) as handle:
            loaded = json.load(handle)
        return _check_threshold_keys(loaded, "filter.thresholds_file")
    inline = filter_section.get("thresholds")
    if isinstance(inline, dict):
        return _check_threshold_keys(inline, "filter.thresholds")
    if required:
        raise ConfigError(
            "No cutoffs yet. Open report.html, move the sliders, download "
            "thresholds.json, and set filter.thresholds_file to that file."
        )
    return None

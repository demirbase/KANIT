// Shared helpers of the KANIT processes. Every step writes into the shared results
// tree (config/config.yaml paths); a task emits receipt.json, which the tasks that
// depend on it stage as input, so that -resume reruns them when it reruns.

// the steps read --config_overlay through KANIT_CONFIG_OVERLAY, set only when there
// is one (an empty env value makes Nextflow warn on every task)
def py(String script) {
    String overlay = params.config_overlay ?
        "KANIT_CONFIG_OVERLAY='${new File(params.config_overlay.toString()).getAbsolutePath()}' " : ''
    "${overlay}${params.python} ${projectDir}/scripts/${script}"
}

// receipt.json of a task, which is also its step manifest (scripts/receipt.py): the
// process, attempt and resources of the task; tools maps a tool's name to the
// command that prints its version
def receipt(task, String step, key, Map tools = [:]) {
    def memory = task.memory ? " --memory '${task.memory}'" : ''
    def versions = tools.collect { name, command -> " --tool '${name}=${command}'" }.join('')
    "${py('receipt.py')} --step ${step} --key '${key}' --process '${task.process}' " +
        "--attempt ${task.attempt} --cpus ${task.cpus}${memory}${versions} > receipt.json"
}

// the receipt of a task that runs on the login node outside the containers
def hostReceipt(task, String step, key) {
    "${params.host_python} ${projectDir}/scripts/receipt.py --step ${step} --key '${key}' " +
        "--process '${task.process}' --attempt ${task.attempt} --cpus ${task.cpus} > receipt.json"
}

def stubReceipt(task, String step, key) {
    "${py('receipt.py')} --step ${step} --key '${key}' --process '${task.process}' --stub " +
        "> receipt.json"
}

def meta(String organism, String antibiotic) {
    [id: "${organism}__${antibiotic}", organism: organism, antibiotic: antibiotic]
}

// config/config.yaml with --config_overlay merged over it, read for the sizes of the
// fan-outs (the steps read the same through KANIT_CONFIG_OVERLAY)
def kanitConfig() {
    def yaml = new org.yaml.snakeyaml.Yaml()
    Map cfg = yaml.load(new File(params.config_yaml).text)
    if (params.config_overlay) {
        cfg = deepMerge(cfg, (yaml.load(new File(params.config_overlay.toString()).text) ?: [:]) as Map)
    }
    cfg
}

Map deepMerge(Map base, Map over) {
    Map out = new LinkedHashMap(base)
    over.each { k, v ->
        out[k] = (v instanceof Map && out[k] instanceof Map) ? deepMerge(out[k] as Map, v as Map) : v
    }
    out
}

// The commands of a batch side by side, each with its share of the task's cores
// (--threads), each writing to its own log. The task fails with the highest exit status
// of its commands; their logs' last lines go to its error output. A retry, after the job
// ran out of memory or time, runs half as many at once as the attempt before, in waves;
// what the commands finished before is kept (--skip-done, the chunk keys).
def parallel(task, List items, Closure command) {
    int wave = Math.max(1, Math.ceil(items.size() / Math.pow(2, task.attempt - 1)) as int)
    def waves = items.withIndex().collate(wave).collect { part ->
        part.collect { item, i -> "${command(item)} --threads \$T > part_${i}.log 2>&1 &\npids+=(\$!)" }
            .join('\n') + '\nwait_all'
    }
    """\
T=\$(( ${task.cpus} / ${wave} )); [ \$T -ge 1 ] || T=1
rc=0
pids=()
wait_all() { for p in "\${pids[@]}"; do c=0; wait \$p || c=\$?; [ \$c -gt \$rc ] && rc=\$c; done; pids=(); }
${waves.join('\n')}
if [ \$rc -ne 0 ]; then tail -n 20 part_*.log >&2; exit \$rc; fi"""
}

// Shared helpers of the KANIT processes. Every step writes into the shared results
// tree (config/config.yaml paths); a task emits receipt.json, which the tasks that
// depend on it stage as input, so that -resume reruns them when it reruns.

def py(String script) {
    "${params.python} ${projectDir}/scripts/${script}"
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

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

def stubReceipt(task, String step, key) {
    "${py('receipt.py')} --step ${step} --key '${key}' --process '${task.process}' --stub " +
        "> receipt.json"
}

def meta(String organism, String antibiotic) {
    [id: "${organism}__${antibiotic}", organism: organism, antibiotic: antibiotic]
}

// config/config.yaml, read for the sizes of the fan-outs (the steps read it themselves)
def kanitConfig() {
    new org.yaml.snakeyaml.Yaml().load(new File(params.config_yaml).text)
}

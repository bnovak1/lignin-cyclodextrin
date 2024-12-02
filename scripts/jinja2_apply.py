from jinja2 import Environment, FileSystemLoader
import sys

# Get the template file name, output format from the command line
template_file = sys.argv[1]
output_format = sys.argv[2]

# Define custom delimiters
custom_env = Environment(
    loader=FileSystemLoader('.'),
    comment_start_string='<!--',
    comment_end_string='-->'
)

# Load the template
template = custom_env.get_template(template_file)

# Define the context (variables) for the template
context = {
    'output_format': output_format
}

# Render the template with the context
rendered_content = template.render(context)

# Save the rendered content to a new file
output_file = f"{output_format}.qmd"
with open(output_file, 'w', encoding='utf-8') as f:
    f.write(rendered_content)
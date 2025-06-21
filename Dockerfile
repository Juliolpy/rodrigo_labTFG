# Use Python 3.12 slim image as base
FROM python:3.12-slim

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Poetry globally
RUN curl -sSL https://install.python-poetry.org | POETRY_HOME=/opt/poetry python3 - && \
    ln -s /opt/poetry/bin/poetry /usr/local/bin/poetry

# Set working directory
WORKDIR /app

# Copy poetry configuration files and README.md first
COPY pyproject.toml poetry.lock* README.md ./

# Copy source code before installing dependencies
COPY src/ ./src/

# Copy the wrapper script
COPY columbo-cli /usr/local/bin/columbo-cli

# Make the wrapper script executable
RUN chmod +x /usr/local/bin/columbo-cli

# Configure Poetry to not create virtual environment and install in system Python
RUN poetry config virtualenvs.create false && \
    poetry config virtualenvs.in-project false

# Install dependencies and the project itself with scripts
RUN poetry install --no-interaction --no-ansi

# Copy data files
COPY data/ ./data/

# Create welcome script
RUN echo '#!/bin/bash\n\
echo "==============================================="\n\
echo "🚀 Welcome to COLUMBO Design Docker Container!"\n\
echo "==============================================="\n\
echo ""\n\
echo "This container contains the COLUMBO diagnostic system"\n\
echo "for designing genetic parts and molecular beacons."\n\
echo ""\n\
echo "📋 Available commands:"\n\
echo "  • columbo-cli --help                    - Show help and usage"\n\
echo "  • columbo-cli <fasta_file>              - Analyze FASTA file (JSON output)"\n\
echo "  • columbo-cli <fasta_file> --output pickle - Analyze FASTA file (Pickle output)"\n\
echo "  • python -m columbo_design.cli --help   - Alternative CLI invocation"\n\
echo "  • poetry run columbo-cli --help         - Run via Poetry"\n\
echo ""\n\
echo "📁 Data files available in /app/data/"\n\
echo "📁 Source code available in /app/src/"\n\
echo ""\n\
echo "🔬 Example usage:"\n\
echo "  columbo-cli /app/data/test.fasta"\n\
echo "  columbo-cli /app/data/mers.fasta --output pickle"\n\
echo ""\n\
echo "==============================================="\n\
echo ""\n\
exec "$@"' > /usr/local/bin/welcome.sh && chmod +x /usr/local/bin/welcome.sh

# Create a non-root user for security
RUN useradd --create-home --shell /bin/bash app && \
    chown -R app:app /app
USER app

# Set welcome script as entrypoint
ENTRYPOINT ["/usr/local/bin/welcome.sh"]

# Default command is bash for interactive access
CMD ["/bin/bash"] 
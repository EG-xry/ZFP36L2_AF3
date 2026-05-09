#!/bin/bash

# AlphaFold-Metainference Docker Helper Script
# This script helps you build and run the Docker container easily

set -e  # Exit on error

cd "$(dirname "$0")"

# Colors for output
GREEN='\033[0.32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}=================================="
echo "AlphaFold-Metainference Docker Helper"
echo -e "==================================${NC}"
echo ""

# Function to check if Docker is running
check_docker() {
    if ! docker info > /dev/null 2>&1; then
        echo -e "${RED}ERROR: Docker is not running!${NC}"
        echo "Please start Docker Desktop and try again."
        exit 1
    fi
}

# Function to build the image
build_image() {
    echo -e "${YELLOW}Building Docker image... This will take 20-30 minutes.${NC}"
    echo "You can monitor progress in another terminal with:"
    echo "  docker-compose logs -f"
    echo ""
    docker-compose build
    echo -e "${GREEN}Build complete!${NC}"
}

# Function to start the container
start_container() {
    echo -e "${YELLOW}Starting container...${NC}"
    docker-compose up -d
    echo -e "${GREEN}Container started!${NC}"
    echo ""
    echo "To enter the container, run:"
    echo "  $0 shell"
}

# Function to stop the container
stop_container() {
    echo -e "${YELLOW}Stopping container...${NC}"
    docker-compose down
    echo -e "${GREEN}Container stopped!${NC}"
}

# Function to enter the container shell
enter_shell() {
    echo -e "${GREEN}Entering container shell...${NC}"
    echo "Type 'exit' to leave the container"
    echo ""
    docker-compose exec afmi-tutorial /workspace/start.sh
}

# Function to show status
show_status() {
    echo "Docker Status:"
    docker-compose ps
    echo ""
    if docker-compose ps | grep -q "Up"; then
        echo -e "${GREEN}✓ Container is running${NC}"
        echo ""
        echo "Available commands:"
        echo "  $0 shell    - Enter the container"
        echo "  $0 stop     - Stop the container"
        echo "  $0 logs     - View container logs"
    else
        echo -e "${YELLOW}Container is not running${NC}"
        echo ""
        echo "Available commands:"
        echo "  $0 start    - Start the container"
        echo "  $0 build    - Rebuild the container"
    fi
}

# Function to show logs
show_logs() {
    docker-compose logs -f
}

# Function to clean up
cleanup() {
    echo -e "${YELLOW}Cleaning up Docker resources...${NC}"
    docker-compose down --rmi all --volumes
    echo -e "${GREEN}Cleanup complete!${NC}"
}

# Main menu
check_docker

case "${1:-}" in
    build)
        build_image
        ;;
    start)
        start_container
        ;;
    stop)
        stop_container
        ;;
    shell|bash|sh)
        enter_shell
        ;;
    restart)
        stop_container
        start_container
        ;;
    status)
        show_status
        ;;
    logs)
        show_logs
        ;;
    cleanup|clean)
        cleanup
        ;;
    *)
        echo "Usage: $0 {build|start|stop|shell|restart|status|logs|cleanup}"
        echo ""
        echo "Commands:"
        echo "  build    - Build the Docker image (first time setup, ~20-30 min)"
        echo "  start    - Start the container"
        echo "  stop     - Stop the container"
        echo "  shell    - Enter the container shell"
        echo "  restart  - Restart the container"
        echo "  status   - Show container status"
        echo "  logs     - Show container logs"
        echo "  cleanup  - Remove all Docker resources"
        echo ""
        echo "Quick start:"
        echo "  1. $0 build     # First time only"
        echo "  2. $0 start     # Start container"
        echo "  3. $0 shell     # Enter container"
        echo ""
        show_status
        exit 1
        ;;
esac


define([
    'jquery',
    'app/ui/layout',
    'app/ui/controls',
    'backbone',
    'view/thread/list',
    'view/thread/create-modal',
    'bootstrap'
], function($, layout, controls, Backbone, ThreadListView, ThreadCreateModal) {

    /**
     * Thread list
     */
    var threadList = new ThreadListView();
    threadList.render();

    /**
     * Thread create modal
     */
    $('#thread-create-modal').on('click', function() {
        var threadCreateModal = new ThreadCreateModal();
        threadCreateModal.render();
    });

});

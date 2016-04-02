define([
    'jquery',
    'backbone',
    'handlebars',
    'collection/thread/thread',
    'text!templates/thread/list.html'
], function ($, Backbone, Handlebars, ThreadCollection, template) {

    var ThreadListView = Backbone.View.extend({
        el: $("#thread-list-container"),

        template: Handlebars.compile(template),

        initialize: function () {
            this.collection = new ThreadCollection();
            this.collection.add({title: "Thread #1 title"});
            this.collection.add({title: "Thread #2 title"});
        },

        render: function() {
            this.$el.html(this.template({threads: this.collection.toJSON()}));

            return this;
        }
    });

    return ThreadListView;
});
